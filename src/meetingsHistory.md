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
    2. Do the subsampling convergence tests for protein mean, pearson'r pval vs gls' pval
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
        4. Why is the Effect Size Mostly Positive, Biologically?
    9. Test another interaction model, where the small model is given by Py ~ Px + M + Drug, only on the model with FDR and no PCA. Check if there is an improvement on the analysability of the volcano Plot
    10. Repeat the analysis in 6 to the volcano Plot with PCA and FDR correction and Compare
    11. (Last Resource if nothing comes out of the previous analyses, because it is more cumbersome) Fazer versão deste modelo de interacção Linear Mixed Models com matrix de covariância da VAE Proteomics
    12. Do the Interaction Model for the Gene Essentiallty instead of Drug Response:
        1.  Get the GeneEsssentially (Pvalues of the respective effect sizes) and Effect sizes of each gene towards every cell line
        2.  Ajustar PValue por Coluna com multiple testing, per Gene in this case
        3.  Convert the ids of Broad to Sanger like The Professor said on Google Chat, use the id mapping provided
        4.  Filtrar genes para que fiquemos com entre 5 a 6k
            1.  Escolher genes que sejam essenciais para pelo menos x cell lines, ou seja, que tenha menos que y p-value e estejam presentes em pelo menos x cell lines (É um guessing game)
    13. Questiona TUDO!
    14. Começa a organizar a escrita, Títulos <-> Sub-títulos <-> Bullet Points
    15. Escreve apenas aquilo que pensas que se manterá imutado até à entrega da tese, datasets e whatnot
    16. FAZ FÉRIAS E PAUSAS
    17. ~~Find a way to correct the fdr interaction model, so that we have a way to descriminate either drugs or PPI's. A work arround could be to use the non inputed data ... choose one (PPI or Drug) and explain the choice. But a drawback from this is that drugs with low samples would get the least penalization, while we are performing the same number of tests in all drugs, maybe the fdr model should be discontinued...~~

    Findings:
    18. The hat is composed of the same PPI in the same direction Py ~ Px, with 686 drugs, where Py is MRPS30 and Px is MRPL41
       1. não faz sentido ser devido ao n (nº de samples em comum entre droga e PPI e respectiva inflação do pValue), pois as associações de sentido inverso também estariam sobre-representadas.
       2. Uma hipotese é número de drogas que podem fazer associação com o PPI (686 drogas), porque este valor afectaria apenas o fdr e possivelmente no modelo inicial onde estamos a fazer a correção na globalidade este fdr cairia acima do threshold, e agora como a penalização é menos signicativa já está abaixo do threshold. Não obstante uma penalização de 686 associações é bastante musculado no quadro global de penalizações cujo percentil 75% é de 221
       3. Outra hipotese é que tem um p-value alto, e um effect-size minusculo devido a estarmos a testar a significância da droga e não só do termo de interação, logo a droga é relevante para explicar Py porque a Drug response aumenta um pouco com aumentos de Py, mas muito provavelmente o termo de interação em si não deve ser muito relevante para explicar a Py, daí o effect size minusculo, ao fazer o modelo novo com Drug no small model este chapéu deve desaparecer.
    19. Some good examples were found on the 7/8, these constitue mainly two complexes:
        1.  PSM (PSMA and PSMB)->This complex plays numerous essential roles within the cell by associating with different regulatory particles;forms the 26S proteasome and thus participates in the ATP-dependent degradation of ubiquitinated protein;plays a key role in the maintenance of protein homeostasis by removing misfolded or damaged proteins that could impair cellular functions, and by removing proteins whose functions are no longer required;mediates ubiquitin-independent protein degradation
        2.  RPS ->
            1. RPSA  -> This gene encodes a high-affinity, non-integrin family, laminin receptor 1. This receptor has been variously called 67 kD laminin receptor, 37 kD laminin receptor precursor (37LRP) and p40 ribosome-associated protein. Laminins, a family of extracellular matrix glycoproteins, are the major noncollagenous constituent of basement membranes. They have been implicated in a wide variety of biological processes including cell adhesion, differentiation, migration, signaling, neurite outgrowth and metastasis. It has been observed that the level of the laminin receptor transcript is higher in colon carcinoma tissue and lung cancer cell line than their normal counterparts. Also, there is a correlation between the upregulation of this polypeptide in cancer cells and their invasive and metastatic phenotype. 
            2. RPS9 -> The protein belongs to the S4P family of ribosomal proteins. It is located in the cytoplasm. Variable expression of this gene in colorectal cancers compared to adjacent normal tissues has been observed, although no correlation between the level of expression and the severity of the disease has been found. 
            3. RPS20 -> This gene encodes a ribosomal protein that is a component of the 40S subunit. The protein belongs to the S10P family of ribosomal proteins. It is located in the cytoplasm. This gene is co-transcribed with the small nucleolar RNA gene U54, which is located in its second intron. 
            4. RPL3  -> This gene encodes a ribosomal protein that is a component of the 60S subunit. The protein belongs to the L3P family of ribosomal proteins and it is located in the cytoplasm. The protein can bind to the HIV-1 TAR mRNA, and it has been suggested that the protein contributes to tat-mediated transactivation. This gene is co-transcribed with several small nucleolar RNA genes, which are located in several of this gene's introns. Alternate transcriptional splice variants, encoding different isoforms, have been characterized. 
            5. RPL7 -> This gene encodes a ribosomal protein that is a component of the 60S subunit. The protein belongs to the L30P family of ribosomal proteins. It contains an N-terminal basic region-leucine zipper (BZIP)-like domain and the RNP consensus submotif RNP2. In vitro the BZIP-like domain mediates homodimerization and stable binding to DNA and RNA, with a preference for 28S rRNA and mRNA. The protein can inhibit cell-free translation of mRNAs, suggesting that it plays a regulatory role in the translation apparatus. It is located in the cytoplasm. The protein has been shown to be an autoantigen in patients with systemic autoimmune diseases, such as systemic lupus erythematosus.
    20. The subsampling tests showed that using the same dataset, Protein Mean 75PV, and using the same proxy value for Protein Interaction, the GLS needed less samples to recall correctly most of the PPI that the full sample model could. With the first model with a subsampling of 5 where the mean of the n models recalled more than 90% of the full model's PPI, the homologous in the Pearson model was below this threshold for the respective full model. In the GLS subsampling run we can see a faster covergence and less variance on the AUC calculated for each of the n models for each subsampling number. Thus, the usefulness of GLS as a method of recalling PPI with fewer samples should be further researched.