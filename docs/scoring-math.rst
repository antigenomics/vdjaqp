Scoring formalism
=================

The notation follows the manuscript and its Supplementary Note 2.
All logarithms in the scoring model are natural.

Error calibration and repertoire background
--------------------------------------------

For true residue :math:`u`, called residue :math:`v`, and error probability :math:`e`, the substitution matrix :math:`S` is conditional on an error. The indicator equals one when its condition holds.

.. math::
   :label: S1

   E_e(u,v)=(1-e)\mathbf1[u=v]+eS_{uv}.

The prior weight :math:`W_{uv}` uses residue frequency :math:`\pi_v`, monoisotopic masses :math:`m_u,m_v`, boost :math:`\beta`, and mass scale :math:`\tau`. Normalize over alternatives :math:`b` to obtain :math:`M`. Combine this prior with substitution counts :math:`N` using pseudocount strength :math:`\alpha`.

.. math::
   :label: S2

   \begin{aligned}
   W_{uv} &= \pi_v\bigl[1+\beta\exp(-|m_u-m_v|/\tau)\bigr],\quad u\ne v, \\
   M_{uv} &= \frac{W_{uv}}{\sum_{b\ne u}W_{ub}},\quad u\ne v,\qquad M_{uu}=0, \\
   S_{uv} &= \frac{N_{uv}+\alpha M_{uv}}{\sum_b N_{ub}+\alpha}.
   \end{aligned}

For cell :math:`h`, :math:`C_{hv}` counts background residues of type :math:`v`, :math:`n_h=\sum_v C_{hv}`, and :math:`\rho_v` is the pooled residue frequency. Cells :math:`j` in backoff group :math:`g(h)` contribute to the stabilized frequency :math:`B_{hv}` with pseudocount strength :math:`\lambda`.

.. math::
   :label: S3

   B_{hv}=\frac{\sum_{j:g(j)=g(h)}C_{jv}+\lambda\rho_v}{\sum_{j:g(j)=g(h)}n_j+\lambda},
   \qquad U_{hv}=C_{hv}+\lambda B_{hv}.

Slot :math:`s` belongs to cell :math:`h(s)` and contains residue :math:`x_s`. Its indicator :math:`k_s` is one if it contributed to the background. Subtract that contribution, apply numerical floor :math:`\epsilon_b`, and normalize by :math:`\nu_s`.

.. math::
   :label: S4

   \begin{aligned}
   \nu_s &= n_{h(s)}-k_s+\lambda, \\
   f_s(v) &= \frac{\max\{U_{h(s)v}-k_s\mathbf1[x_s=v],\epsilon_b\}}{\nu_s}, \\
   g_s(v) &= \frac{\max\{\sum_u U_{h(s)u}S_{uv}-k_sS_{x_sv},\epsilon_b\}}{\nu_s}.
   \end{aligned}

Placement and residue evidence
------------------------------

For placement :math:`a`, sum over valid aligned positions :math:`V_a`. Each slot :math:`s_i` has candidate residue :math:`x_{s_i}` and the peptide reports :math:`y_{pi}` with calibrated error :math:`e_{pi}`. The candidate call model is :math:`E_e`. The background mixes the context-specific residue frequency :math:`f_s` and its error-transformed distribution :math:`g_s`.

.. math::
   :label: S5

   \begin{aligned}
   \mathcal B_s(v;e) &= (1-e)f_s(v)+e\,g_s(v), \\
   \ell_a &= \sum_{i\in V_a}\log
   \frac{E_{e_{pi}}(x_{s_i},y_{pi})}{\mathcal B_{s_i}(y_{pi};e_{pi})}.
   \end{aligned}

For observation :math:`o`, all placements in :math:`A_o` compete with a unit-weight unassigned option. For peptide :math:`p`, the tilde quantities are bounded confidence, rank and proposal count, respectively. The resulting responsibility is :math:`\gamma_a`.

.. math::
   :label: S6

   w_a=e^{\ell_a}\frac{\tilde z_p^2}{\tilde b_p\tilde d_p},
   \qquad \gamma_a=\frac{w_a}{1+\sum_{b\in A_o}w_b}.

At slot :math:`s=(c,r)`, the agreeing placements :math:`A_s^+` use I/L equivalence. The numerical floor is :math:`\epsilon_u`; :math:`I_{\max}` caps residue information in nats.

.. math::
   :label: S7

   Q_s=1-\prod_{a\in A_s^+}\max(1-\gamma_a,\epsilon_u),
   \qquad I_s=\min[-\log f_s(x_s),I_{\max}].

The junction scoring interval of candidate :math:`c` is :math:`\mathcal J_c`. The negative placements are :math:`\mathcal A_c^-` and :math:`D_{\min}` is the contradiction floor. The total, junction, contradiction and combined evidence are :math:`H_c,J_c,D_c,F_c`, respectively.

.. math::
   :label: S8

   \begin{aligned}
   H_c &= \sum_r Q_{cr}I_{cr}, &
   J_c &= \sum_{r\in\mathcal J_c} Q_{cr}I_{cr}, \\
   D_c &= \max\!\left(D_{\min},\sum_{a\in\mathcal A_c^-}\gamma_a\ell_a\right), &
   F_c &= H_c+D_c.
   \end{aligned}

When both isotype groups occur, the final score uses prior weights :math:`p_G,p_O` and candidate counts :math:`n_G,n_O` for the IgG group :math:`G` and the other group :math:`O`. With one group, :math:`T_c=F_c`.

.. math::
   :label: S9

   T_c=F_c+\mathbf1[c\in G]\log\frac{p_G/n_G}{p_O/n_O}.

Redundancy analysis
-------------------

For a peptide call :math:`p` of length :math:`L_p`, let :math:`\mathcal M_{pk}` be valid placements in clonotype :math:`k`. The conditional call probability :math:`P(y_{pi}\mid x_{ai})` uses the fitted error model at position :math:`i` in placement :math:`a`. Empty placement sets contribute zero to the mean across all :math:`K` clonotypes.

.. math::
   :label: S10

   L_{pk}=\max_{a\in\mathcal M_{pk}}\prod_{i=1}^{L_p} P(y_{pi}\mid x_{ai}),
   \qquad \overline{L}_p=\frac{1}{K}\sum_{k=1}^{K}L_{pk}.

For spectrum :math:`s`, :math:`\eta_{sk}` is the largest log ratio to the common peptide likelihood across the alternative calls, with the confidence, rank and proposal-count factor. :math:`\omega_s` is the inverse number of spectra with the same I/L-equivalent alternative set. Fit the nonnegative weights :math:`\theta_k`, including the unassigned weight :math:`\theta_0`, on the simplex.

.. math::
   :label: S11

   \mathcal Q(\theta)=\sum_s \omega_s\log\!\left(\theta_0+\sum_{k=1}^{K}\theta_k\exp(\eta_{sk})\right),
   \qquad \sum_{k=0}^{K}\theta_k=1.

The loss compares full and reduced fitted objectives :math:`\widehat{\mathcal Q}_F,\widehat{\mathcal Q}_{-k}`. Subtract the remaining improvement bound :math:`\Delta_{-k}=\max(0,\max_j g_j-\sum_s\omega_s)`, where :math:`g_j` is the reduced-objective gradient with respect to remaining weight :math:`\theta_j`.

.. math::
   :label: S12

   B_k=\max(0,\widehat{\mathcal Q}_F-\widehat{\mathcal Q}_{-k}-\Delta_{-k}).

For the separate manuscript redundancy analysis, :math:`B_c` measures the
loss of fit after removing a candidate, and :math:`\kappa(n_c)` is reliability
at RNA support :math:`n_c`. A positive removal score supports retention.

.. math::
   :label: removal-score

   S_c = \max\{0, B_c + \log\kappa(n_c)\}.
