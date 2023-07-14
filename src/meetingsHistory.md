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
    4. Remove the 3 dots in the top 2 of the model, that have a high value for x and low for y and drugRes = 1, redo the linear model, see if the interaction effect size decreases (Those points were creating a relationship where there was none)
    5. Criar um PPI set com Corum == 1 e fdr < 1%, mandar nº para o Emanuel. Tens de utilizar os arrays dos pValues com multipletesting e voilá, fazer o mesmo para o pearson pvals correlation matrices não vá o diabo tecê-las
    6. De forma ao scatter plot ser mais palpável, e para perceber melhor as respostas extremas de certas cell lines a drogas, fazer um joint kernel plot? Ou um scatter com com as cores em contínuo, ou seja a resposta à droga não binarizada
    7. De forma a combater biases que estão relacionados com o número de vezes que o evento ocorre (nº de samples onde está presente aquelas duas proteínas), fazer multiple testing por ppi no interaction model, por exemplo as proteínas do proteosoma são bastante abundantes e é um complexo com imensos elemntos o que faz com que as suas correlações sejam muito pouco específicas e desinteressantes (Perguntar melhor que biases são estes em detalhe)
    8. Outra forma de combater esses biases é regress out as 10 PC do PCA mais fortes na esperança de explicarem os confouding factors que estão a criar um interaction term onde ele não existe, colocariamos no modelo da interacção em vez das growth_props como . Fazer barplot de dos PC's (a sua variância explicada) e o cumulativo desta variãncia.
    9. (Outlook) Fazer versão deste modelo de interacção Linear Mixed Models com matrix de covariância da VAE Proteomics
    
    Carry Over:
       1. Recalculate the AUC for the various pairwise models
       2. Do the subsampling convergence tests for protein mean, pearson'r pval vs gls' pval