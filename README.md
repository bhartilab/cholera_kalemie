# Minimal data and code for the analyses and figures in "Impact of a multi-pronged cholera intervention in an endemic setting"

# Structure

```markdown
├── R code
│   ├── mcmc_sampler_migr.R
│   ├── mcmc_sampler_migrS.R
│   ├── mcmc_sampler_migrN.R
│   ├── other_functions_migr.R
│   ├── other_functions_migrS.R
│   ├── other_functions_migrN.R
│   └── figures.R
│
└── data
    ├── cases.rda
    ├── cum_df.rda
    ├── cum_df2.rda
    ├── fit_df.rda
    ├── num_env.rda
    ├── num_est.rda
    ├── prop.rda
    ├── foi_env.rda
    ├── growth_env.rda
    ├── env_data.rda
    ├── vacci_res_c1.rda
    ├── vacci_step3_s1.rda
    ├── vacci_step3_s2.rda
    ├── vacci_step3_s3.rda
    └── vacci_step3_s4.rda
```

# R code

Below is an overview of the scripts and their functions. Additional details are provided in notes in the code. All the packages necessary are either loaded at the beginning of the script or specified along the code.

We provide the code to sample from the posterior distribution for

* the full model with both susceptible and infected individuals migrating (*mcmc_sampler_migr.R*)

* the model with only susceptible indivuals migrating (*mcmc_sampler_migrS.R*)

* the model with no migration (*mcmc_sampler_migrN.R*)

The code for the figures uses the samples and estimates from the last one: the model with no migration. It is the model presented in the publication.

## *mcmc_sampler_migr.R*

This script includes the code use for sampling from the posterior distribution for the model assuming that both susceptible and infected individuals migrate. You need to choose between the set of assumptions we explored (see the Supplementary Information (SI)) on:

* vaccine effectiveness (see more in '*vacci_step3_s1* to *vacci_step3_s4*' in *data* below)

* the penalty term applied to account for the spatially targeted nature of the vaccination ('*s_capt*')

* the ratio between susceptible and infectious among the mobile population ('*ratio_S_I*')

The sampler requires loading other functions from the script *other_functions_migr.R* :

* *dXdt*

* *log_post_theta.1*

* *log_post_theta.2*

* *log_post_theta.3*

Running the sampler requires a vector *cases* including the weekly aggregated number of suspected cholera cases. These data can be made available by contacting Klaudia Porten (Klaudia.PORTEN@epicentre.msf.org).

You must define the total number of samples (*iterations*), the burn in period (*burn.in*), the batch size to save the chain regularly (*write.states*).

The values for the environmental variables are loaded at the beginning of the script. More information on them in the SI. The preprocessed variables are available as rda files (see below in *data*).

## *mcmc_sampler_migrS.R*

This script includes the code use for sampling from the posterior distribution for the model assuming that only susceptible individuals migrate. You need to choose between the set of assumptions we explored (see the Supplementary Information (SI)) on:

* vaccine effectiveness (see more in '*vacci_step3_s1* to *vacci_step3_s4*' in *data* below)

* the penalty term applied to account for the spatially targeted nature of the vaccination ('*s_capt*')

The sampler requires loading other functions from the script *other_functions_migrS.R* :

* *dXdt*

* *log_post_theta.1*

* *log_post_theta.2*

* *log_post_theta.3*

Running the sampler requires a vector *cases* including the weekly aggregated number of suspected cholera cases. These data can be made available by contacting Klaudia Porten (Klaudia.PORTEN@epicentre.msf.org).

You must define the total number of samples (*iterations*), the burn in period (*burn.in*), the batch size to save the chain regularly (*write.states*).

The values for the environmental variables are loaded at the beginning of the script. More information on them in the SI. The preprocessed variables are available as rda files (see below in *data*).

## *mcmc_sampler_migrN.R*

This script includes the code use for sampling from the posterior distribution for the model assuming that there is no migration. You need to choose between the set of assumptions we explored (see the Supplementary Information (SI)) on:

* vaccine effectiveness (see more in '*vacci_step3_s1* to *vacci_step3_s4*' in *data* below)

* the penalty term applied to account for the spatially targeted nature of the vaccination ('*s_capt*')

The sampler requires loading other functions from the script *other_functions_migrN.R* :

* *dXdt*

* *log_post_theta.1*

* *log_post_theta.2*

* *log_post_theta.3*

Running the sampler requires a vector *cases* including the weekly aggregated number of suspected cholera cases. These data can be made available by contacting Klaudia Porten (Klaudia.PORTEN@epicentre.msf.org).

You must define the total number of samples (*iterations*), the burn in period (*burn.in*), the batch size to save the chain regularly (*write.states*).

The values for the environmental variables are loaded at the beginning of the script. More information on them in the SI. The preprocessed variables are available as rda files (see below in *data*).

## *other_functions_migr.R*

This script includes functions necessary to run the sampler and additional functions producing the estimates to assess the impact of the intervention and explore the relative contributions of the environmental reservoir and seasonal migration. It is used by *mcmc_sampler_migr.R*.

Those functions rely on the sample from the posterior distribution saved after ring the sampler. We provide preprocessed statistics allowing a reproduction of all our results (see in *data* below and in the code) and the code used to get them (provided as notes). Those functions include:

* the functions used to pre-process the environmental variables: *rad_lag2*, *rain_lag2*, *sst_lag*, and *chlor_lag* (more information in the SI). Those functions process external data freely available. We provide information on where to download them in the code, but we do not provide them because we do not own them.

* a function to reconstruct the chain from the files saved after running the sampler: *build_sample2*

* a function to produce the estimates to assess the impact of the intervention: *fit_impact_intervention*. It creates a file *estimate.rda* used to create pre-processed estimates necessary to produce the figures. Those pre-processed estimates are provide in */data*.

* a function to produce the estimates to explore alternative vaccinatio strategies: *alternative_vaccination_strategies*. It creates a file *vacci_res_c1.rda*, provided in */data*.

* a function to produce estimates to assess the relative contributions of the environmental reservoir and seasonal mobility: *migration_environment*. It creates a file *contrib.rda* used to create pre-processed estimates necessary to produce the figures. Those pre-processed estimates are provide in */data*. 

## *other_functions_migrS.R*

This script includes functions necessary to run the sampler and additional functions producing the estimates to assess the impact of the intervention and explore the relative contributions of the environmental reservoir and seasonal migration. It is used by *mcmc_sampler_migrS.R*.

Those functions rely on the sample from the posterior distribution saved after ring the sampler. We provide preprocessed statistics allowing a reproduction of all our results (see in *data* below and in the code) and the code used to get them (provided as notes). Those functions include:

* the functions used to pre-process the environmental variables: *rad_lag2*, *rain_lag2*, *sst_lag*, and *chlor_lag* (more information in the SI). Those functions process external data freely available. We provide information on where to download them in the code, but we do not provide them because we do not own them.

* a function to reconstruct the chain from the files saved after running the sampler: *build_sample2*

* a function to produce the estimates to assess the impact of the intervention: *fit_impact_intervention*. It creates a file *estimate.rda* used to create pre-processed estimates necessary to produce the figures. Those pre-processed estimates are provide in */data*.

* a function to produce the estimates to explore alternative vaccinatio strategies: *alternative_vaccination_strategies*. It creates a file *vacci_res_c1.rda*, provided in */data*.

* a function to produce estimates to assess the relative contributions of the environmental reservoir and seasonal mobility: *migration_environment*. It creates a file *contrib.rda* used to create pre-processed estimates necessary to produce the figures. Those pre-processed estimates are provide in */data*. 

## *other_functions_migrN.R*

This script includes functions necessary to run the sampler and additional functions producing the estimates to assess the impact of the intervention and explore the relative contributions of the environmental reservoir and seasonal migration. It is used by *mcmc_sampler_migrN.R*.

Those functions rely on the sample from the posterior distribution saved after ring the sampler. We provide preprocessed statistics allowing a reproduction of all our results (see in *data* below and in the code) and the code used to get them (provided as notes). Those functions include:

* the functions used to pre-process the environmental variables: *rad_lag2*, *rain_lag2*, *sst_lag*, and *chlor_lag* (more information in the SI). Those functions process external data freely available. We provide information on where to download them in the code, but we do not provide them because we do not own them.

* a function to reconstruct the chain from the files saved after running the sampler: *build_sample2*

* a function to produce the estimates to assess the impact of the intervention: *fit_impact_intervention*. It creates a file *estimate.rda* used to create pre-processed estimates necessary to produce the figures. Those pre-processed estimates are provide in */data*.

* a function to produce the estimates to explore alternative vaccinatio strategies: *alternative_vaccination_strategies*. It creates a file *vacci_res_c1.rda*, provided in */data*.

* a function to produce estimates to assess the relative contributions of the environmental reservoir and seasonal mobility: *migration_environment*. It creates a file *contrib.rda* used to create pre-processed estimates necessary to produce the figures. Those pre-processed estimates are provide in */data*. 

## *figures.R*

This script includes the code required to produce the figures in the manuscript. They are provided as functions, *fig1*, *fig2*, *fig3*, and *fig4*. They use pre-processed data sets provided in */data*.

Those functions produce eps files that we then impoved on adobe illustrator. This second step is not covered here.

# Data

The folder */data* includes minimal data sets necessary to reproduce all the figues of the manuscript (Figures 1 to 4).

They include interim results after processing the posterior distribution of the parameters to visually assess the model fit and its estimates (*fit_df.rda*), estimate the impact of the intervention, alternative scenarios (*prop.rda*, *num_est.rda*, *cum_df.rda*, and *vacci_res_c1.rda*), and the contribution of the seasonal migration and the environmental reservoir (*foi_env.rda*, *growth_env.rda*, *env_data.rda*, *num_env.rda*, and *cum_df2.rda*). They also include minimal data sets used to explore alternative assumptions on vaccine effectiveness (*vacci_step3_s1.rda* to *vacci_step3_s4.rda*).

The code used to produce those interim results is available in *R code* (see below).

When external data sets are used, we provide information on where to find them, but we do not provide them because we are not the original owners.

## *cases*

It is a single vector of length 118 (number of weeks) with the number of suspected cholera cases in the only cholera treatment center of Kalemie aggregated by week.

It is used in Figure 2 and by the mcmc sampler (see in the section dedicated to the R code below).

This file is **NOT** freely available because they are sureveillance data. It can be made available by contacting Klaudia Porten (Klaudia.PORTEN@epicentre.msf.org).

## *vacci_step3_s1* to *vacci_step3_s4*

Those files include vectors of length 118 used to define the step function *eta* used in the mcmc sampler while using different assumptions on the vaccine effectiveness.

More on the assumptions is available in the Supplementary Information (SI) in sections *Step function $\eta$*, *Vaccine coverage*, and *Vaccine effectiveness*.

*vacci_step3_s1* to *vacci_step3_s4* correspond to VE assumptions 1 to 4 in Table S1 of the SI.

## Pre-processed estimates to assess the fit visually and explore its predictions of the evolution of the compartments of the model

It only includes *fit_df.rda*.

It is necessary to produce Figure 2 and is produced by functions available in */R code/other_functions_migrN.R*.

## Pre-processed estimates to assess the impact of the intervention and alternative vaccination strategies

Those pre-processed estimates include:

* *prop.rda*

* *num_est.rda*

* *cum_df.rda*

* *vacci_res_c1.rda*

They are necessary to produce Figure 3 and are produced by functions available in */R code/other_functions_migrN.R*.


## Pre-processed estimates to assess the relative contributions of the environmental reservoir and seasonal migration

Those pre-processed estimates include:

* *foi_env.rda*

* *growth_env.rda*

* *env_data.rda*

* *num_env.rda*

* *cum_df2.rda*

They are necessary to produce Figure 4 and are produced by functions available in */R code/other_functions_migrN.R*.
