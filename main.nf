nextflow.enable.dsl=2

process ANALYSE {
    tag "$sample"
    publishDir { "${params.outdir}/${sample}" }, mode: 'copy', overwrite: false
    input:
    tuple val(sample), path(candidates), path(peptides), val(producer), path(model)
    output:
    tuple val(sample), path('analysis')
    script:
    """
    '${params.python}' '${projectDir}/scripts/workflow_analysis.py' '${candidates}' '${peptides}' '${model}' analysis --producer '${producer}' --cpus ${task.cpus}
    """
}

process SIMULATE_AIRR {
    publishDir params.outdir, mode: 'copy', overwrite: false
    input:
    path repertoire
    path panel
    output:
    path 'simulation'
    path 'simulation-resources.json'
    script:
    def rssLimit = params.mode == 'joint-airr' ? 8000000000 : 15000000000
    """
    '${params.python}' '${projectDir}/scripts/profile_build.py' --output "\$PWD/simulation-resources.json" --project-root '${projectDir}/..' --workload-root "\$PWD/simulation" --timeout 1650 --rss-limit ${rssLimit} --disk-limit 10000000000 --cpus ${task.cpus} -- '${params.python}' '${projectDir}/scripts/workflow_simulation.py' --repertoire "\$PWD/${repertoire}" --target-ids "\$PWD/${panel}" --source-format '${params.producer}' --output "\$PWD/simulation" --physical-scope '${params.physical_scope}' --ms2-per-bin ${params.ms2_per_bin} --cpus ${task.cpus}
    """
}

process SYNTHETIC {
    publishDir params.outdir, mode: 'copy', overwrite: false
    input:
    path references
    output:
    path 'synthetic'
    script:
    """
    PYTHONPATH='${params.arda_overlay}' '${params.python}' '${projectDir}/scripts/smoke_all.py' --output synthetic --reference-dir '${references}' --cpus ${task.cpus}
    """
}

process REAL_COHORT {
    publishDir params.outdir, mode: 'copy', overwrite: false
    input:
    path input_root
    output:
    path 'real-analysis'
    script:
    """
    '${params.python}' '${projectDir}/scripts/workflow_real.py' --input-root '${input_root}' --output real-analysis --cpus ${task.cpus}
    """
}

process POPULATION {
    publishDir params.outdir, mode: 'copy', overwrite: false
    output:
    path 'population'
    path 'population-resources.json'
    script:
    def baseArg = params.references ? "--reference-dir '${file(params.references).toAbsolutePath()}'" : ''
    def rssLimit = params.mode == 'joint' ? 7000000000 : 15000000000
    def referenceArg = params.additional_references ? "--additional-reference-dir '${file(params.additional_references).toAbsolutePath()}'" : ''
    """
    PYTHONPATH='${params.arda_overlay}' '${params.python}' '${projectDir}/scripts/profile_build.py' --output "\$PWD/population-resources.json" --project-root '${projectDir}/..' --workload-root "\$PWD/population" --timeout 1650 --rss-limit ${rssLimit} --disk-limit 10000000000 --cpus ${task.cpus} -- '${params.python}' '${projectDir}/scripts/workflow_population.py' --output "\$PWD/population" --cpus ${task.cpus} --total-antibodies ${params.antibodies} --memory-fraction ${params.memory_fraction} --pairs-per-chain ${params.pairs_per_chain} --mutation-rate ${params.mutation_rate} --protein-concentration ${params.protein_concentration} --seed ${params.seed} --read-seed ${params.read_seed} --gradient-seconds ${params.gradient_seconds} --ms2-per-bin ${params.ms2_per_bin} --purpose resource-profile --stage-timeout ${params.stage_timeout} ${referenceArg} ${baseArg}
    """
}

// Paths cross a shell boundary. Reject shell metacharacters before submission.
def safe(value) {
    if (value == null || !(value.toString() ==~ /[A-Za-z0-9_., \/:+=-]+/)) error "Unsafe or missing workflow argument"
    return value
}

workflow {
    [projectDir, params.python, params.outdir, params.producer].each { safe(it) }
    if (file(params.outdir).exists()) error 'Output directory exists; choose a fresh --outdir'
    if (!(params.physical_scope in ['all','panel'])) error 'Invalid physical scope'
    if (!(params.ms2_per_bin.toString() ==~ /[1-9][0-9]*/)) error 'ms2_per_bin must be positive'
    if (params.mode in ['simulate-airr','joint-airr'] && !(params.ms2_per_bin.toString() in ['5','10'])) error 'Supplied-AIRR cases support ms2_per_bin 5 or 10'
    if (!(params.cpus.toString() ==~ /[1-9][0-9]*/)) error 'cpus must be positive'
    if (workflow.resume) error 'Resume is disabled: fresh calculations are required.'
    if (params.analysis_executor != null) {
        if (params.mode != 'analysis') error 'Analysis resource profiles support --mode analysis only'
        if (!(params.analysis_task_cpus.toString() ==~ /[1-9][0-9]*/)) error 'analysis_task_cpus must be positive'
        if ((params.analysis_task_memory as nextflow.util.MemoryUnit).toBytes() >
            (params.analysis_memory_budget as nextflow.util.MemoryUnit).toBytes())
            error 'analysis_memory_budget must cover one analysis task'
    }
    if (params.mode == 'analysis') {
        if (!params.samples) error 'Provide --samples samples.csv'
        def seen = new HashSet()
        rows = Channel.fromPath(params.samples, checkIfExists: true).splitCsv(header: true).map { r ->
            if (!(r.sample ==~ /[A-Za-z0-9_-]+/) || !(r.producer ==~ /[A-Za-z0-9_. -]+/)) error 'Invalid sample or producer'
            if (!seen.add(r.sample)) error 'Duplicate sample name'
            [r.candidates, r.peptides, r.model].each { safe(it) }
            tuple(r.sample, file(r.candidates, checkIfExists:true), file(r.peptides, checkIfExists:true), r.producer, file(r.model, checkIfExists:true))
        }
        ANALYSE(rows.ifEmpty { error 'Empty sample sheet' })
    } else if (params.mode in ['simulate-airr','joint-airr']) {
        [params.repertoire, params.panel].each { safe(it) }
        if (params.mode == 'joint-airr') {
            def allocation=Math.min(params.cpus as int,(System.getenv('SLURM_CPUS_PER_TASK') ?: Runtime.runtime.availableProcessors()) as int)
            if (allocation < 8) error 'joint-airr requires at least 8 allocated cores'
            safe(params.input_root)
            REAL_COHORT(file(params.input_root,checkIfExists:true))
        }
        SIMULATE_AIRR(file(params.repertoire, checkIfExists:true), file(params.panel, checkIfExists:true))
    } else if (params.mode == 'synthetic') {
        safe(params.arda_overlay)
        if (!file(params.arda_overlay,checkIfExists:true).isDirectory()) error 'Run sh envs/arda/setup_overlay.sh first'
        safe(params.references)
        SYNTHETIC(file(params.references, checkIfExists:true))
    } else if (params.mode == 'real-cohort') {
        safe(params.input_root)
        REAL_COHORT(file(params.input_root, checkIfExists:true))
    } else if (params.mode in ['population','joint']) {
        safe(params.arda_overlay)
        if (!file(params.arda_overlay, checkIfExists:true).isDirectory()) error 'Run sh envs/arda/setup_overlay.sh first'
        if (params.references) {
            safe(params.references)
            if (!file(params.references,checkIfExists:true).isDirectory()) error 'references must be a raw public-reference directory'
        }
        if (params.mode == 'joint') {
            def allocation=Math.min(params.cpus as int,(System.getenv('SLURM_CPUS_PER_TASK') ?: Runtime.runtime.availableProcessors()) as int)
            if (allocation < 8) error 'joint requires at least 8 allocated cores'
            safe(params.input_root)
            REAL_COHORT(file(params.input_root,checkIfExists:true))
        }
        if (params.additional_references) {
            safe(params.additional_references)
            if (!file(params.additional_references, checkIfExists:true).isDirectory()) error 'additional_references must be a raw public-reference directory'
        }
        ['antibodies','pairs_per_chain'].each { name ->
            if (!(params[name].toString() ==~ /[1-9][0-9]*/)) error "${name} must be a positive integer"
        }
        if (!(params.seed.toString() ==~ /[0-9]+/) || !(params.read_seed.toString() ==~ /[0-9]+/)) error 'seeds must be nonnegative integers'
        ['memory_fraction','mutation_rate','gradient_seconds','stage_timeout'].each { name ->
            def value=params[name].toString()
            if (!(value ==~ /[0-9]+(?:[.][0-9]+)?(?:[eE][+-]?[0-9]+)?/) || !Double.isFinite(value.toDouble())) error "${name} must be finite numeric input"
        }
        if (!(params.protein_concentration in ['equal','head80'])) error 'Invalid protein_concentration'
        // Existing Python CLI validates biological ranges/family divisibility.

        POPULATION()
    } else error 'mode must be analysis, simulate-airr, synthetic, real-cohort, population, joint or joint-airr'
}
