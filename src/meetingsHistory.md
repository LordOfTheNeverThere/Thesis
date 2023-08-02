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
    5.  ~~Understand why there is a hat on the volcano plot of FDR correction without PCA~~ -- The hat is composed of the same PPI in the same direction Py ~ Px, with all the drugs, where Py is MRPS30 and Px is MRPL41, não faz sentido se dever nem ao n (nº de samples em comum entre droga e PPI e respectiva inflação do pValue), pois as associações de sentido inverso também estariam sobre-representadas ou número de drogas que podem fazer associação com o PPI (686 drogas), porque este valor afectaria apenas o fdr e não o p-value do log likelihood test.
    6.  ~~Maybe have in the scatterplots the drug response given by Size and not Colour~~
    7.  Correlacionar resíduos do small model Py ~ Px + M, com Drogas
    8.  Analyse in concrete the Volcano Plot of th FDR without PCA (Remember the Professor said that we need only one example of biologival interest in order to make our model viable)
        1.  Select certain points or clusters of points
        2.  Look into the Fudge Factor, those weird tagents with seconf order polynomial shape. As a method to select points of interest
        3.  Colorir Volcano com várias cores: (O objectivo é fugir de de linhas polinomiais contínuos, ou seja uma variável não deve estar associada a uma projecção específica do volcano plot, deve ser all over the place, combinar isto talvez com fudge factor)
            1.  Colorir com #Samples
            2.  Px ou Py
            3.  Drug
            4.  Etc
        4. Why is the Effect Size Mostly Positive, Biologically?
    9. Repeat the analysis in 6 to the volcano Plot with PCA and FDR correction and Compare
    10. Test another interaction model, where the small model is given by Py ~ Px + M + Drug, only on the model with FDR and no PCA. Check if there is any difference between that and the homologous, by doing the same analysis done in 8 and 9
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

            