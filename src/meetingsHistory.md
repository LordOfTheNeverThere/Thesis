5/7/23
Tasks:
    1. Zscore the X variable in all regressions - DONE
    2. Always have intercept on all linear methods, otherwise we might be forcing a line whihc does not explain all the variability in data in the best way - Done
    3. In a linear model check if there is Normality and Homoscedasticity, (?independence of samples?), Optional
    4. Redo the dead residuals model -  DONE (no major improvements, not adding growth_props increases R^2 to 0.11 with Lung and Heamo)
    5. Redo the TLS regression class, solve the problem with the -Vxy/Vzz -DONE
    6. Focus then on the interaction term Linear Model -DONE
        1. Save the number of common samples, n - DONE
        2. Save the effectsize of the interaction term, how much drug response alters the explanatory power of Px -DONE
        3. Do Log ratio test to see if the interaction term and Drug Response term are significant py ~ M + px + drug + px*drug -DONE
        4. Do a volcano Plot and then a scatter plot of the samples with highest interaction term in abs and highest Y axis coordinate
    7. Recalculate the AUC for the various pairwise models
    8. Do the subsampling convergence tests for protein mean, pearson'r pval vs gls' pval


Doubts:
    1. Should I add random normal noise to binary confoundind factors like we are doing in ResidualsLinearModel? -Done
    2. In the intercept Linear Model, should I have a regression for each drug and PPI and not account for all the drugs simultaneously? The meaning of the wilk's test is not very meaningful -Done


12/7/23
Tasks:
    1. ~~Send the 60 head of the results of the interaction linear model~~
    2. ~~Normalise the drug in the interaction linear model~~
    3. ~~Remove the black border in the alrealdy translucent dots of the volcano plot~~
    4. ~~Remove the 3 dots in the top 2 of the model, that have a high value for x and low for y and drugRes = 1, redo the linear model, see if the interaction effect size decreases (Those points were creating a relationship where there was none, this could be do to outliers which had a very high drugResponse, which were also outliers in the high drug response set)~~
    5. ~~Criar um PPI set com Corum == 1 e fdr < 1%, mandar nº para o Emanuel. Tens de utilizar os arrays dos pValues com multipletesting e voilá, fazer o mesmo para o pearson pvals correlation matrices não vá o diabo tecê-las~~
    6. ~~De forma ao scatter plot ser mais palpável, e para perceber melhor as respostas extremas de certas cell lines a drogas, fazer um joint kernel plot? Ou um scatter com com as cores em contínuo, ou seja a resposta à droga não binarizada, ver se o Emanuel Gosta, talvez adicionar para ali uns boxplots idk!~~
    7. ~~De forma a combater biases que estão relacionados com o número de vezes que o evento ocorre (nº de samples onde está presente aquelas duas proteínas), fazer multiple testing por ppi no interaction model, por exemplo as proteínas do proteosoma são bastante abundantes e é um complexo com imensos elemntos o que faz com que as suas correlações sejam muito pouco específicas e desinteressantes (this could be do to outliers which had a very high drugResponse, which were also outliers in the high drug response set)~~
    8. ~~Outra forma de combater esses biases é regress out as 10 PC do PCA (vaeProteomics) mais fortes na esperança de explicarem os confouding factors que estão a criar um interaction term onde ele não existe, colocariamos no modelo da interacção em vez das growth_props como . Fazer barplot de dos PC's (a sua variância explicada) e o cumulativo desta variãncia.~~
    9. ~~(Outlook) Fazer versão deste modelo de interacção Linear Mixed Models com matrix de covariância da VAE Proteomics~~
    10. ~~See if Emanuel Likes the volcano plot as it is along with the scatterPlot!!!~~
    11. ~~Update: Recalculate the AUC of the various Pairwise Models~~
        1.  ~~Redo the barplots~~
        2.  ~~Do an additional one, with the a grid like structure and multiple subplots one for each Dataset, inquery Emanuel on more details on this~~
    
    Carry Over:
       1. ~~Do the subsampling convergence tests for protein mean, pearson'r pval vs gls' pval~~
    Doubts:tenta
        2. (In the interaction model) Why is the X not only the interaction term instead of both interaction term and, a resposta do stackoverflow que baseou esta tentativa refere apenas em testar o termo de interação https://stats.stackexchange.com/questions/435644/is-there-a-method-to-look-for-significant-difference-between-two-linear-regressi  This makes sense since we should test if the linear model is better explained by the interaction term than without it. So we want to see a model where the interaction between drugResponse and Px is statistically important at explaining Py. With our current model we might have low pValues because the Drug Response acompanies Py, therefore there is no interaction present. Thus we are not associating the alteration of PxPy correlation with a change in Drug Response 


26/7/23
Tasks:
    1. ~~Do Roc Curve of VAE vs Original in terms of pearson'r (PRIORITY)~~
    2. ~~Do the subsampling convergence tests for protein mean, pearson'r pval vs gls' pval~~
    3. ~~Update: Recalculate the AUC of the various Pairwise Models~~
        1.  ~~Redo the barplots~~
        2.  ~~Do an additional one, with the a grid like structure and multiple subplots one for each Dataset, inquery Emanuel on more details on this~~
    4. ~~The Bar Plot in the Scree Plot is Correct but the cummulative is of the explained_variance and not explained_variance_ratio, so we want the variance explained by each PC for that data, regardless of of the variance explained by all the PC's~~
    5.  ~~Understand why there is a hat on the volcano plot of FDR correction without PCA~~
    6.  ~~Maybe have in the scatterplots the drug response given by Size and not Colour~~
    7.  Correlacionar resíduos do small model Py ~ Px + M, com Drogas (Method has been made but in order for it to be used the Interaction Models need to be recalculated)
    8.  Analyse in concrete the Volcano Plot of th FDR without PCA (Remember the Professor said that we need only one example of biological interest in order to make our model viable)
        1.  Select certain points or clusters of points
        2.  Look into the Fudge Factor, those weird tagents with seconf order polynomial shape. As a method to select points of interest
        3.  ~~Colorir Volcano com várias cores: (O objectivo é fugir de de linhas polinomiais contínuos, ou seja uma variável não deve estar associada a uma projecção específica do volcano plot, deve ser all over the place, combinar isto talvez com fudge factor)~~
            1.  ~~Colorir com #Samples~~
            2.  ~~Px ou Py~~
            3.  ~~Drug~~
            4.  ~~Etc~~
        4. ~~Why is the Effect Size Mostly Positive, Biologically?~~
    9. ~~Test another interaction model, where the small model is given by Py ~ Px + M + Drug, only on the model with FDR and no PCA. Check if there is an improvement on the analysability of the volcano Plot~~
    10. Repeat the analysis in 8 to the volcano Plot with PCA and FDR correction and Compare
    11. (Last Resource if nothing comes out of the previous analyses, because it is more cumbersome) Fazer versão deste modelo de interacção Linear Mixed Models com matrix de covariância da VAE Proteomics
    12. (Optional) Test another test statistic to check if the interaction term is significant
    13. Do the Interaction Model for the Gene Essentiallty instead of Drug Response:
        1.  Get the GeneEsssentially (Pvalues of the respective effect sizes) and Effect sizes of each gene towards every cell line
        2.  Ajustar PValue por Coluna com multiple testing, per Gene in this case
        3.  Convert the ids of Broad to Sanger like The Professor said on Google Chat, use the id mapping provided
        4.  Filtrar genes para que fiquemos com entre 5 a 6k
            1.  Escolher genes que sejam essenciais para pelo menos x cell lines, ou seja, que tenha menos que y p-value e estejam presentes em pelo menos x cell lines (É um guessing game)
    14. ~~Why are the large residues of both models different? Shouldn't they be the reflection of the same model, maybe something is going on the classes.py~~ I should calculate ssr instead of sum of the resiudals
    15. Questiona TUDO!
    16. Começa a organizar a escrita, Títulos <-> Sub-títulos <-> Bullet Points
    17. Escreve apenas aquilo que pensas que se manterá imutado até à entrega da tese, datasets e whatnot
    18. FAZ FÉRIAS E PAUSAS
    19. ~~Find a way to correct the fdr interaction model, so that we have a way to descriminate either drugs or PPI's. A work arround could be to use the non inputed data ... choose one (PPI or Drug) and explain the choice. But a drawback from this is that drugs with low samples would get the least penalization, while we are performing the same number of tests in all drugs, maybe the fdr model should be discontinued...~~

    Findings:
    20. The hat is composed of the same PPI in the same direction Py ~ Px, with 686 drugs, where Py is MRPS30 and Px is MRPL41
       1. não faz sentido ser devido ao n (nº de samples em comum entre droga e PPI e respectiva inflação do pValue), pois as associações de sentido inverso também estariam sobre-representadas.
       2. Uma hipotese é número de drogas que podem fazer associação com o PPI (686 drogas), porque este valor afectaria apenas o fdr e possivelmente no modelo inicial onde estamos a fazer a correção na globalidade este fdr cairia acima do threshold, e agora como a penalização é menos signicativa já está abaixo do threshold. Não obstante uma penalização de 686 associações é bastante musculado no quadro global de penalizações cujo percentil 75% é de 221
       3. Outra hipotese é que tem um p-value alto, e um effect-size minusculo devido a estarmos a testar a significância da droga e não só do termo de interação, logo a droga é relevante para explicar Py porque a Drug response aumenta um pouco com aumentos de Py, mas muito provavelmente o termo de interação em si não deve ser muito relevante para explicar a Py, daí o effect size minusculo, ao fazer o modelo novo com Drug no small model este chapéu deve desaparecer.
    21. Some good examples were found on the 7/8, these constitue mainly two complexes:
        1.  PSM (PSMA and PSMB)->This complex plays numerous essential roles within the cell by associating with different regulatory particles;forms the 26S proteasome and thus participates in the ATP-dependent degradation of ubiquitinated protein;plays a key role in the maintenance of protein homeostasis by removing misfolded or damaged proteins that could impair cellular functions, and by removing proteins whose functions are no longer required;mediates ubiquitin-independent protein degradation
        2.  RPS ->
            1. RPSA  -> This gene encodes a high-affinity, non-integrin family, laminin receptor 1. This receptor has been variously called 67 kD laminin receptor, 37 kD laminin receptor precursor (37LRP) and p40 ribosome-associated protein. Laminins, a family of extracellular matrix glycoproteins, are the major noncollagenous constituent of basement membranes. They have been implicated in a wide variety of biological processes including cell adhesion, differentiation, migration, signaling, neurite outgrowth and metastasis. It has been observed that the level of the laminin receptor transcript is higher in colon carcinoma tissue and lung cancer cell line than their normal counterparts. Also, there is a correlation between the upregulation of this polypeptide in cancer cells and their invasive and metastatic phenotype. 
            2. RPS9 -> The protein belongs to the S4P family of ribosomal proteins. It is located in the cytoplasm. Variable expression of this gene in colorectal cancers compared to adjacent normal tissues has been observed, although no correlation between the level of expression and the severity of the disease has been found. 
            3. RPS20 -> This gene encodes a ribosomal protein that is a component of the 40S subunit. The protein belongs to the S10P family of ribosomal proteins. It is located in the cytoplasm. This gene is co-transcribed with the small nucleolar RNA gene U54, which is located in its second intron. 
            4. RPL3  -> This gene encodes a ribosomal protein that is a component of the 60S subunit. The protein belongs to the L3P family of ribosomal proteins and it is located in the cytoplasm. The protein can bind to the HIV-1 TAR mRNA, and it has been suggested that the protein contributes to tat-mediated transactivation. This gene is co-transcribed with several small nucleolar RNA genes, which are located in several of this gene's introns. Alternate transcriptional splice variants, encoding different isoforms, have been characterized. 
            5. RPL7 -> This gene encodes a ribosomal protein that is a component of the 60S subunit. The protein belongs to the L30P family of ribosomal proteins. It contains an N-terminal basic region-leucine zipper (BZIP)-like domain and the RNP consensus submotif RNP2. In vitro the BZIP-like domain mediates homodimerization and stable binding to DNA and RNA, with a preference for 28S rRNA and mRNA. The protein can inhibit cell-free translation of mRNAs, suggesting that it plays a regulatory role in the translation apparatus. It is located in the cytoplasm. The protein has been shown to be an autoantigen in patients with systemic autoimmune diseases, such as systemic lupus erythematosus.
    22. The subsampling tests showed that using the same dataset, Protein Mean 75PV, and using the same proxy value for Protein Interaction, the GLS needed less samples to recall correctly most of the PPI that the full sample model could. With the first model with a subsampling of 5 where the mean of the n models recalled more than 90% of the full model's PPI, the homologous in the Pearson model was below this threshold for the respective full model. In the GLS subsampling run we can see a faster covergence and less variance on the AUC calculated for each of the n models for each subsampling number. Thus, the usefulness of GLS as a method of recalling PPI with fewer samples should be further researched.
    23. Why is the interaction effect size mostly positive: *Show a powerpoint as a way to show this information in the best way possible*
        1.  If the effect size is positive it means that the interaction between drug res and Px amplifies the impact on the explanatory power towards Py. It indicates that the combined effect is greater than the sum of their individual contributions, so the Px affects the value of Drug Res and vice versa
        2.  On the other hand, if the effect size is negative the interaction term mitigates the individual effects of each of the predictor vars, so it acts as a buffering mechanism... The increase of Px decreases the explanatory power of DrugRes towards Py and vice versa
        3.  Since the effect size is mostly positive, it means two things:
            1.  Py ~ Px * (beta_Px + drugRes*beta_int) + drugRes*beta_drugRes, The increase in Drug Response magnifies the impact of Px in perdicting Py. So the fact is that the Drug Response helps Px explain Py better. Thence, the correlation between Py and Px increases, a complex is more likely formed.
            2.  Py ~ Px * beta_Px + drugRes*(beta_drugRes + Px*beta_int),The increase in Px magnifies the explanatory power of drugRes towards Py. Thence, the correlation between drugRes and Py increases when we increase Px, this is not really interesting. Since if we look at the equation above knowing that Px and Py come from a set of known ppis where the majority are positively correlated, so the effect size of the interaction term is sure to be positive since most Px are correlated positively with Py, and thus add a little explanatory power to that of the drugRes.
            3.  So in sum only the first effect is actually interesting on a biological sense, since the  drug response impacts the slope between Px and Py and thus complex formation
    24. In terms of memory the output data from the interaction model occupies about 313 bytes per row, so in the gene dependency model if we have 5k genes we will have an outputing dataframe of about numGenes*numPPIs*2 = 5000*12726*2 = 127.260.000 rows... which would cost 39884102592.77743 bytes or 38Gb, if we instead use 2k we get a much more modest size of 15Gb
    25. While converting the Broad Ids to Sanger Ids using the mapping file the prof provided, these Broad Ids had no Sanger Matches... Is it normal? 54 samples in total which 5% of total samples
        1.  IndexError: ACH-000467 not found in mappingFile
            IndexError: ACH-000833 not found in mappingFile
            IndexError: ACH-001063 not found in mappingFile
            IndexError: ACH-001172 not found in mappingFile
            IndexError: ACH-001370 not found in mappingFile
            IndexError: ACH-001393 not found in mappingFile
            IndexError: ACH-001481 not found in mappingFile
            IndexError: ACH-001543 not found in mappingFile
            IndexError: ACH-001675 not found in mappingFile
            IndexError: ACH-001691 not found in mappingFile
            IndexError: ACH-001719 not found in mappingFile
            IndexError: ACH-001834 not found in mappingFile
            IndexError: ACH-001839 not found in mappingFile
            IndexError: ACH-001844 not found in mappingFile
            IndexError: ACH-001862 not found in mappingFile
            IndexError: ACH-001961 not found in mappingFile
            IndexError: ACH-001986 not found in mappingFile
            IndexError: ACH-001990 not found in mappingFile
            IndexError: ACH-002014 not found in mappingFile
            IndexError: ACH-002024 not found in mappingFile
            IndexError: ACH-002035 not found in mappingFile
            IndexError: ACH-002070 not found in mappingFile
            IndexError: ACH-002084 not found in mappingFile
            IndexError: ACH-002239 not found in mappingFile
            IndexError: ACH-002446 not found in mappingFile
            IndexError: ACH-002458 not found in mappingFile
            IndexError: ACH-002459 not found in mappingFile
            IndexError: ACH-002460 not found in mappingFile
            IndexError: ACH-002461 not found in mappingFile
            IndexError: ACH-002462 not found in mappingFile
            IndexError: ACH-002463 not found in mappingFile
            IndexError: ACH-002464 not found in mappingFile
            IndexError: ACH-002465 not found in mappingFile
            IndexError: ACH-002467 not found in mappingFile
            IndexError: ACH-002471 not found in mappingFile
            IndexError: ACH-002485 not found in mappingFile
            IndexError: ACH-002486 not found in mappingFile
            IndexError: ACH-002508 not found in mappingFile
            IndexError: ACH-002510 not found in mappingFile
            IndexError: ACH-002511 not found in mappingFile
            IndexError: ACH-002512 not found in mappingFile
            IndexError: ACH-002523 not found in mappingFile
            IndexError: ACH-002526 not found in mappingFile
            IndexError: ACH-002531 not found in mappingFile
            IndexError: ACH-002650 not found in mappingFile
            IndexError: ACH-002654 not found in mappingFile
            IndexError: ACH-002693 not found in mappingFile
            IndexError: ACH-002710 not found in mappingFile
            IndexError: ACH-002785 not found in mappingFile
            IndexError: ACH-002799 not found in mappingFile
            IndexError: ACH-002800 not found in mappingFile
            IndexError: ACH-002834 not found in mappingFile
            IndexError: ACH-002847 not found in mappingFile
            IndexError: ACH-002922 not found in mappingFile
            IndexError: ACH-002926 not found in mappingFile


24/8/23
    1. ~~Regarding the hat on the Model that includes the drug on the Large Model we ought to remove it from fruther analyses if it ever arises again~~
    2. ~~Only work with the Large Model Model, that is Py ~ Px + M + Drug/Gene + Drug/Gene*Px, the volcano plot of the Small Model Model is strangely worse... (Well let's see)~~
    3. ~~Do the Synthetic data examples that will help understand the DRInteraction Model:~~
        1. Starting with Positevely correlated Proteins: beta_Px > 0
          1. One where the low IC50 is the linear line and the high IC50 are cluster above -> Expected: beta_int > 0
          2. One where the low IC50 is the linear line and the high IC50 are cluster bellow -> Expected: beta_int < 0
          3. One where the high IC50 is the linear line and the low IC50 are cluster above -> Expected: beta_int < 0
          4. One where the high IC50 is the linear line and the low IC50 are cluster bellow -> Expected: beta_int > 0
        2. Repeat the same with Negatively Correlated Proteins: beta_Px < 0
            One where the low IC50 is the linear line and the high IC50 are cluster above -> Expected: beta_int > 0 
            One where the low IC50 is the linear line and the high IC50 are cluster bellow -> Expected: beta_int < 0
            One where the high IC50 is the linear line and the low IC50 are cluster above -> Expected: beta_int < 0
            One where the high IC50 is the linear line and the low IC50 are cluster bellow -> Expected: beta_int > 0
        3. Write a table with +/- on the betas values, and betas from intercept, Px, Drug/gene, Interaction
    4. ~~Redo the same Drug Interaction Model but now instead of using ppi that are in corum under a specific fdr in the proteomics matrix, use those ppi seen in string900 and biogrid, since corum is ribossome protein enriched which gives us a not so cool bias, since ribossomal complexes contain many proteins. But it is unlikely that they have a big role on cancer response to drugs.~~
       1. Find a way to select just a few if there is a too large number.
    5. ~~Re-Select the trully essential genes, since it seems that the p-value and effect size is not a suficient metric since if we consider genes it atleast one sample < 0.01 (Professor set as most restrictive filtration) we will get most of the genes, so not a reliable way to filter~~
       1. ~~Do a sort of minMax Scalling with the median of the NE and E (NonEssential and Essential, respectively) from the files the professor sent to google chat~~
          1. ~~So per sample we calculate the median of essential genes medianEss and the counterRespective medianNonEss, we then standardise for each value x in column ((x - medianNonEss) / (medianNonEss - medianEss))~~
             1. So 1std will be the difference in previousSTD / (medianNonEss - medianEss)
             2. So essential genes will have standardiseX < 0 since the greater the essentiality the less the log of fold change, since the fold change becomes a less and lesser ratio, less cells survived after testing. Remenber that log fold change is log(Final/Initial)
             3. And non essential genes constitue the postive numbers
             4. If x = medianEss the transformedX = -1
             5. transformedX is only 1 when the x =(2*medianNonEss - medianEss)
          2. ~~Then select genes that have at least one sample that is less than -0.5, or any other value so that we have only genes which are heterogenous in essentiality in at least one sample~~
       2. ~~After the first filtration we used the filtrated genes on the second filtration, using Fisher's Skewdness test~~
          1. Aceepting only genes with skewness in distribution less than -2
          2. We want genes wich are negatively skewed somewhat since that means that they have an handful of observayion where the log fold change is more negative than the majority of the remaining samples. So it is a gene that is non essential for most samples, but for a set of few is essential. These samples are of biological interest since we ought to understand what happens in those samples that makes the cell lines suscepitable to the loss of expression of that gene
       3. The set of genes of this two step filtration should be of interest, since these metric are not purely imperical but have some biological and statistical sense into it, counter respective.
    6. ~~Make scatterplot have first points of smallest Drug Response or Gene Dependency, by sorting as ascending the drug Response.~~


Findings:
    1. By filtering the gene set in the gene dependency problem, If after the median scalling per sample ((x - medianNonEss) / (medianNonEss - medianEss)), making all samples comparable since they were scalled using the same principle. Where we select only genes that have at least one sample less than -0.5. After wards if we select genes that have a Pearson-Fisher skew less than -1.25 we get a set of 976 genes. geneDependency.filterGenes(skewThresh=-1.25)
    Finnished scaling samples per median of essential and non essential set of genes and selecting only genes with at least one sample of value less than -0.5,

    From all the gene data of shape 17931, from the first filtration the set of genes is of size 12360 

    From all the gene data of shape 12360, from the second filtration the set of genes is of size 976
    2. The Hat and all the other hats that were seen in the volcano plot, changing y value and ppi were caused by linear regression instability, the effect size of the growthProps were either on the magnitude of the other features, about e-1 or e-2, or due to instability e+12 or e+13. So I replaced the growthProps, that were also problematic on the previous TLS Linear model, due to the same instable nature, with 5 pc's from PCA, that in sum explained 60% of data variation. Perhaps would be unwise to continue to use Growth props in linear regression, maybe it is its sparsity that causes this variability

13/09/23
    1. ~~In all linear regressions done in order to avoid any numerical instability, make sure that no regressions are made, with features that have no variability, std = 0~~
    2. Justify why the choice of -1.25 as a skewness threshold, do like an histogram and say we selected x = -1.25 or whatever then drawing a line.
    3. ~~Redo the interaction linear model:~~
       1. The best model should be DR ~ M + Px + Py + Px*Py
          1. With this model we can find cases where:
             1. The Drug Response alters the Pearson correlation, which was the aim of the first model thus replacing it in functionality. By testing the interaction term PxPy
             2. If Px alters DR alone, so as Px increases so does DR, this biological mechanism does not correspond to the change of the stoichometry of the Px-Py complex since this test captures changes even if the slope between Px and Py remains the same. So this corresponds to another mechanism
             3. If Py alters Dr alone, idem idem idem
          2. And it will require half the time since we do not need to reverse the PPI's
    4. Regarding the convergence of subsampling of Proeteomics75PV using GLS (Calculating AUC), pick one of the subsampling levels where the whole distribution is above the orange line, and repeat the subsampling and AUC calculation 100 times or more to understand if the distribution closes up to the orange line or not... Since I cannot show that chart and expect for people not to question why that subsampling level has more recalling power than all samples... If this test converged to the orange line we would be safe from the prying eyes of reviewers
    5. ~~Make boxplor/histo of the p-values of the 3 tests in the new model. In order to show what mechanism is mostly associated with the drug response. It seams that PxPy has less small pValues than Px and Py alone~~


29/9/23
 1. WRITEEEEEEE! 
    1. Revise what was written before, and see what must be rewritten and what must be added! The algorithmic part is the one which will need more looking up too
    2. Send each finnished chapter to the Professor, to ease revision
    3. You can't put everything on the thesis, focus on the interaction (The other two terms were alrealdy done extensively) term imagery and meaning, but the leftovers in the annex
       1. An example:
          1. Introduce Model V for DR, saying that you used Corum ppis with a certain threshold
          2. Show volcano plot (The drug one is interesting, put the others in annex)
          3. Say that you did exactly the same for another set of protein, String900 union Biogrid, and fds < 0.01
          4. Show volcano plot (The drug one is interesting, put the others in annex)
          5. Show list of most significant association, like top 10, in both cases, compare them, they should have some similarities, but not completely equal.
          6. Repeat the same process but now for Gene Dependency.   

 2. If model IV will be included in the thesis, understand better if The Professor wants it included!
    1. See if the top 30 associations seen in the model IV are included in model's V associations, this should be the case, since model V is an improvement from V. Report this
 3. ~~Revise the Signs in the pictorial Examples~~
 4. ~~Add a biological meaning to it, by associating the sign of the beta of interaction term or combination of signs of other betas with change in Drug Response (Ideally create an association which is simple and holds up in all cases in the same manner)~~
 5. We have 8 possibilities in our models, see how many association we have for each of the eigth classes (after fdr cut off), then plot a volcano plot with this classes as hue, also do an heatmap
    1. In addition we can add a method that adds a new column in each association that says what class it belongs (Begin writting first, at leat a bit)
 6. ~~The 3 distributions of the three pValues of Model V has an artefact, there is suposed to be an uniform distribution, with an inflation on the small p-values. Professor said to look up the internet to find what is it? (After starting writting)~~
    1. After research you found out that something is wrong with your statistical test, maybe the p-value is being calculated on the wrong distro, revise the code, and then rerun the code, in this step you will need aproxametely 3 days, for the 3 models, seen in http://varianceexplained.org/statistics/interpreting-pvalue-histogram/

6/10/23:
1. Regarding the subsampling artefact, in order to try to find out if it happens when we subsample a medium number of samples due to sampling we might not be including outliers that in the full proteomics would decrease or increase the correlation in absolute value, and making the matter worse the inputting by the mean would mean that we would need more outliershy samples to break this general effect perserved in the mean, do the following:
   1. Get two sets of samples, one with only samples with the lowest missing Values, other with highest missing value rate, within the selected few proteins of the inputted dataset.
   2. For each select two distant subsampling numbers, that share the characteristic of having higher auc than all the matrix of inputted proteomics

20/10/23
1.Do the top10 scatter plots, Py, Px, Drug Response, pick one of the figures for a main figure
2. Understan if the gene essentiality seen with SOX1 is tissue type specific
3. Do some schematics to explain some workflow, like the way the two sub sampling tests were made.
4. Extended abstract é 80% copy paste, não necessita mais do que 1 dia at best
5. Enviar 1st draft by monday