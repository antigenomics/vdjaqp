"""Fixed 255-coordinate AIRR abundance reconstruction definitions; no ABC engine."""
import math
import numpy as np

REGIONS=('FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4')
BLOCKS=(('solution_mixture','trypsin'),('gel_mixture','trypsin'),('solution_mixture','chymotrypsin'))
SIGMA=math.log(10)


def relative_amounts(theta):
    """Eleven orthonormal centered log contrasts -> twelve positive amounts, sum 12."""
    theta=np.asarray(theta,dtype=float)
    if theta.shape!=(11,) or not np.isfinite(theta).all():raise ValueError('eleven finite log contrasts required')
    basis=np.zeros((12,11))
    for k in range(1,12):
        basis[:k,k-1]=1/math.sqrt(k*(k+1));basis[k,k-1]=-k/math.sqrt(k*(k+1))
    z=basis@theta;w=np.exp(z-z.max());w=12*w/w.sum()
    if not np.isfinite(w).all() or np.any(w<=0):raise ValueError('unrepresentable positive weights')
    return w


def target_vector(targets):
    """Use one shared chymotrypsin block; exact private ref1..ref12 ordering."""
    rows=targets.to_dicts() if hasattr(targets,'to_dicts') else targets
    ys=[];coverage=[];denominators=[]
    def one(**keys):
        matches=[r for r in rows if all(r.get(k)==v for k,v in keys.items())]
        if len(matches)!=1:raise ValueError('target key missing or duplicated')
        return matches[0]
    for condition,enzyme in BLOCKS:
        y=one(condition=condition,enzyme=enzyme,reference='all',quantity='compatible_called_species')
        if y['units']!='species' or y['value']!=y['numerator'] or int(y['value'])!=y['value'] or y['value']<0:
            raise ValueError('invalid species target')
        ys.append(y['value']);c=[];d=[]
        for i in range(1,13):
            cr=[];dr=[]
            for region in REGIONS:
                r=one(condition=condition,enzyme=enzyme,reference=f'ref{i}',quantity='residue_union_coverage',region=region)
                n,den=r['numerator'],r['denominator']
                if r['units']!='fraction' or den is None or den<=0 or not 0<=n<=den or not math.isclose(r['value'],n/den):
                    raise ValueError('invalid reference coverage target')
                cr.append(r['value']);dr.append(den)
            c.append(cr);d.append(dr)
        coverage.append(c);denominators.append(d)
    # Duplicate display of the shared experimental acquisition must agree exactly.
    shared={ (r['reference'],r['quantity'],r['region']):(r['value'],r['numerator'],r['denominator'])
            for r in rows if r['condition']=='solution_mixture' and r['enzyme']=='chymotrypsin'}
    gel={ (r['reference'],r['quantity'],r['region']):(r['value'],r['numerator'],r['denominator'])
            for r in rows if r['condition']=='gel_mixture' and r['enzyme']=='chymotrypsin'}
    if shared!=gel:raise ValueError('shared chymotrypsin target aliases disagree')
    result=dict(yield_counts=ys,coverage=coverage,denominators=denominators)
    validate_vector(result)
    return result


def validate_vector(summary):
    y=np.asarray(summary['yield_counts'],dtype=float);c=np.asarray(summary['coverage'],dtype=float);d=np.asarray(summary['denominators'],dtype=float)
    if (y.shape!=(3,) or c.shape!=(3,12,7) or d.shape!=c.shape or not np.isfinite(y).all()
            or not np.isfinite(c).all() or not np.isfinite(d).all() or np.any(y<0) or np.any(y!=np.floor(y))
            or np.any(c<0) or np.any(c>1) or np.any(d<=0) or np.any(d!=np.floor(d))):
        raise ValueError('invalid 255-coordinate summary or denominator')
    return y,c,d


def distance(simulated,observed):
    sy,sc,sd=validate_vector(simulated);oy,oc,od=validate_vector(observed)
    if not np.array_equal(sd,od):raise ValueError('reference region denominators differ')
    return float(np.mean(.5*np.abs(np.log1p(sy)-np.log1p(oy))/math.log(2)+.5*np.mean(np.abs(sc-oc),axis=(1,2))))


class BudgetExhausted(RuntimeError):
    """A hard attempted-draw or aggregate-time bound, including incomplete populations."""


def run_abc(forward,observed,output,*,max_draws=255,wall_seconds=21600.,particles=32,populations=3,seed=314159):
    """Published pyABC engine with a hard model-call guard; never resume history.

    Forward receives twelve relative amounts and a draw index. This callback is
    the only scientific execution boundary; no ABC algorithm is implemented here.
    """
    import json,time
    from pathlib import Path
    from datetime import timedelta
    from importlib.metadata import version
    if version('pyabc')!='0.13.1':raise ValueError('qualified pyABC 0.13.1 required')
    import pyabc
    validate_vector(observed)
    if max_draws<=0 or not math.isfinite(wall_seconds) or wall_seconds<=0:raise ValueError('positive inference budgets required')
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    np.random.seed(seed)
    started=time.monotonic();attempted=0
    def model(params):
        nonlocal attempted
        if attempted>=max_draws or time.monotonic()-started>=wall_seconds:
            raise BudgetExhausted('attempted-draw or aggregate-wall budget exhausted')
        attempted+=1
        theta=[params[f'z{k}'] for k in range(11)]
        amounts=relative_amounts(theta)
        record=dict(draw=attempted,log_contrasts=theta,relative_amounts=amounts.tolist(),status='running')
        path=output/f'attempt-{attempted:04d}.json';path.write_text(json.dumps(record,indent=2)+'\n')
        try:
            result=forward(amounts,attempted,max(0.,wall_seconds-(time.monotonic()-started)))
            record.update(status='complete',distance=distance(result,observed),summary=result)
            return result
        except BaseException as error:
            record.update(status='failed',error_type=type(error).__name__,error=str(error));raise
        finally:path.write_text(json.dumps(record,indent=2)+'\n')
    class StrictTransition(pyabc.MultivariateNormalTransition):
        def fit(self,X,w):
            if len(X)<2:raise ValueError('too few particles for unregularized covariance')
            raw=np.atleast_2d(np.cov(X.values,aweights=w,rowvar=False))
            if raw.shape!=(11,11) or not np.isfinite(raw).all() or np.linalg.eigvalsh(raw).min()<=0:
                raise ValueError('singular/nonfinite unregularized log-contrast covariance')
            super().fit(X,w)
            if not np.isfinite(self.cov).all() or np.linalg.eigvalsh(self.cov).min()<=0:
                raise ValueError('singular/nonfinite log-contrast transition covariance')
    class FlooredMedian(pyabc.MedianEpsilon):
        def __call__(self,t):
            return max(.10,float(super().__call__(t)))
    prior=pyabc.Distribution(**{f'z{k}':pyabc.RV('norm',0,SIGMA) for k in range(11)})
    abc=pyabc.ABCSMC(model,prior,distance_function=distance,population_size=particles,
        sampler=pyabc.sampler.SingleCoreSampler(check_max_eval=True),eps=FlooredMedian(initial_epsilon='from_sample',median_multiplier=1.,weighted=True),
        transitions=StrictTransition(scaling=2.,bandwidth_selector=lambda n,d:1.))
    history=abc.new('sqlite:///'+str((output/'abc.sqlite').resolve()),observed)
    status='failed';error_record=None
    try:
        history=abc.run(minimum_epsilon=.10,max_nr_populations=populations,max_total_nr_simulations=max_draws,
                        max_walltime=timedelta(seconds=wall_seconds))
        last_epsilon=float(history.get_all_populations().query('t >= 0')['epsilon'].iloc[-1]) if history.max_t>=0 else None
        status='complete' if history.max_t+1==populations or (last_epsilon is not None and last_epsilon<=.10) else 'budget_limited'
    except BudgetExhausted as error:
        status='budget_limited';error_record=str(error)
    except BaseException as error:
        error_record=f'{type(error).__name__}: {error}';raise
    finally:
        import csv
        population_ess={}
        for population in range(max(0,history.max_t+1)):
            frame,weights=history.get_distribution(m=0,t=population)
            population_ess[str(population)]=float(1/np.sum(np.square(weights)))
            with (output/f'population-{population}-particles.tsv').open('w') as stream:
                fields=['particle','weight']+[f'z{k}' for k in range(11)]+[f'relative_amount_ref{i+1}' for i in range(12)]
                writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
                for index,(_,row) in enumerate(frame.iterrows()):
                    theta=[float(row[f'z{k}']) for k in range(11)];amounts=relative_amounts(theta)
                    writer.writerow(dict(particle=index,weight=float(weights[index]),
                        **{f'z{k}':theta[k] for k in range(11)},
                        **{f'relative_amount_ref{i+1}':float(amounts[i]) for i in range(12)}))
        report=dict(status=status,attempted_draws=attempted,population_effective_sample_size=population_ess,allowed_draws=max_draws,
            wall_seconds=time.monotonic()-started,allowed_wall_seconds=wall_seconds,
            completed_populations=max(0,history.max_t+1),particles_per_population=particles,
            pyabc_version=version('pyabc'),sigma=SIGMA,epsilon_schedule='weighted MedianEpsilon from initial prior sample, floor0.10',
            epsilon_history={str(int(row['t'])):float(row['epsilon']) for _,row in history.get_all_populations().iterrows() if row['t']>=0},
            initial_prior_calibration_is_counted=True,recovery_used=False,
            first_engineering_draw_is_not_prior_particle=True,error=error_record)
        (output/'inference_report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
