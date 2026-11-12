# Reference library

Papers used to build [../experiment/GLOSSARY.md](../experiment/GLOSSARY.md). The tag in square brackets (e.g. [R01]) is how the glossary cites each paper.

- `.pdf` files are the full paper.
- `.xml` files are the full text in JATS XML from PubMed Central. They open in any browser or text editor.
- `_ABSTRACT_ONLY` means the publisher did not allow a full-text download. Use the link to read it.

| Tag | File | Citation | What we take from it |
|---|---|---|---|
| R01 | 01_Battista2019_... .pdf | Battista, Pearcy, Strickland (2019). *Modeling the prescription opioid epidemic.* Bull. Math. Biol. 81. arXiv:1711.03658 | Compartment structure (S, P, A, R), relapse split into spontaneous (σ) and peer-driven (νRA), baseline yearly rates (Table 1), and the result that a prescription route removes the drug-free equilibrium |
| R02 | 02_WhiteComiskey2007_... .pdf | White & Comiskey (2007). *Heroin epidemics, treatment and ODE modelling.* Math. Biosci. 208 | First SIR-type drug model. Initiation by contact with users, treatment, and relapse driven by contact with users |
| R04 | 04_Ruhm2018_... .pdf | Ruhm (2018). *Deaths of despair or drug problems?* NBER WP 24188 | The drug environment (availability and price) explains more of the county death trend than economic decline does. Justifies the dynamic availability layer |
| R05 | 05_Hollingsworth2017_... .pdf | Hollingsworth, Ruhm, Simon (2017). *Macroeconomic conditions and opioid abuse.* NBER WP 23192 / J. Health Econ. | +1 percentage point of county unemployment gives +3.6% opioid death rate. Used for elasticity matching |
| R06 | 06_Monnat2018_... .xml | Monnat (2018). *Factors associated with county-level differences in U.S. drug-related mortality rates.* AJPM 54(5) | Economic and family distress raise mortality; religious establishments are protective; no average rural/urban gap |
| R07 | 07_Monnat2019_... ABSTRACT_ONLY | Monnat et al. (2019). *Using census data to understand county-level differences in overall drug mortality and opioid-related mortality by opioid type.* AJPH. [PMC6611117](https://pmc.ncbi.nlm.nih.gov/articles/PMC6611117/) | Census predictors differ by opioid type (prescription vs heroin vs synthetic) |
| R08 | 08_Pitt2018_... ABSTRACT_ONLY | Pitt, Humphreys, Brandeau (2018). *Modeling health benefits and harms of public policy responses to the US opioid epidemic.* AJPH. [PMC6137764](https://pmc.ncbi.nlm.nih.gov/articles/PMC6137764/) | Compartmental policy model; supply restriction pushes people toward illicit opioids |
| R09 | 09_Ciccarone2019_... .xml | Ciccarone (2019). *The triple wave epidemic: supply and demand drivers of the US opioid overdose crisis.* IJDP | Three waves (prescription, heroin, fentanyl), each with its own supply and demand drivers. Basis for the two-channel availability (Rx + illicit) and the potency variable |
| R10 | 10_Ciccarone2017_... .xml | Ciccarone (2017). *US heroin in transition* (editorial). IJDP | Fentanyl as a supply shock that changes the risk environment, not user behavior |
| R11 | 11_Cerda2022_... ABSTRACT_ONLY | Cerdá et al. (2022). *A systematic review of simulation models to track and address the opioid crisis.* Epidemiol. Rev. [PMC9005056](https://pmc.ncbi.nlm.nih.gov/articles/PMC9005056/) | Map of existing models. Agent-based models are only 16% of them, and almost none are spatial over real density |
| R12 | 12_Sordo2017_... .xml | Sordo et al. (2017). *Mortality risk during and after opioid substitution treatment.* BMJ 357:j1550 | Overdose mortality of 12.7 per 1,000 person-years outside methadone treatment and 2.6 inside it. Baseline for μ_H |
| R13 | 13_Mednick2010_... .pdf | Mednick, Christakis, Fowler (2010). *The spread of sleep loss influences drug use in adolescent social networks.* PLoS ONE 5(3) | A friend's drug use raises own use by about 110%. Evidence for a contact term |
| R14 | 14_AkterDeardon2024_... .pdf | Akter & Deardon (2024). *Conditional logistic individual-level models of spatial infectious disease dynamics.* arXiv:2409.02353 | Spatial individual-level model (ILM) with a distance kernel: P(infect) = 1 − exp(−α Σ d^−β) |
| R15 | 15_RahulDeardon2024_... .pdf | Rahul & Deardon (2024). *Nonparametric spatial kernels for ILMs.* arXiv:2405.00835 | Step-function (piecewise) kernels. Our fixed-radius kernel is the one-step special case |
| R16 | 16_Queen2022_... .pdf | Queen et al. (2022). *Agent-based dynamics of a SPAHR opioid model on social network structures.* arXiv:2202.12261 | Agent-based version of the Battista model. Network structure changes the outcome compared with mean-field ODEs |
| R18 | 18_Tiwari2024_... .pdf | Tiwari et al. (2024). *Measuring network dynamics of opioid overdose deaths in the US.* arXiv:2410.17496 | Deaths in socially connected counties predict own deaths. 97% of people in a withdrawal program knew the person they initiated with |
| R19 | 19_Zibbell2023_... .xml | Zibbell et al. (2023). *Associations between opioid overdose deaths and drugs confiscated by law enforcement and submitted to crime laboratories.* | Each NFLIS fentanyl exhibit per 100k gives +2.97% overdose deaths per quarter. Proxy for illicit potency |

These were cited but could not be downloaded (blocked):

- Bilden, Kumar, Roberts (2025). *Evaluation and replication of a social contagion component of OUD using agent-based simulation.* medRxiv. [link](https://www.medrxiv.org/content/10.1101/2025.02.03.25321377v1). An agent-based model in which proximity to users in the household and neighborhood reproduces the exponential rise in deaths. **Closest prior work to our design.**
- Ruhm (2018). *Drivers of the fatal drug epidemic* (draft). Same conclusion as R04, at county level.
- Rhodes (2009). *Risk environments and drug harms: a social science for harm reduction approach.* IJDP 20(3). Conceptual basis for "environment × individual" rate modifiers.
- Bass (1969). *A new product growth model for consumer durables.* Management Science 15(5). The diffusion equation we reuse for fentanyl penetration.
