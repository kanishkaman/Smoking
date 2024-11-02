import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import f
 
# Loading and preprocessing the data
def load_and_preprocess_data(file_path):
    data = pd.read_csv(file_path, delimiter='\t')
   
    expression_data = np.exp(data.iloc[:, 1:49].values) # Taking the exponent, because the data is log-transformed.
    # expression_data = np.exp(data.iloc[:, 1:49].values)     --> USE THIS, IF YOU WANT TO RUN DIRECTLY ON THE DATA PROVIDED.
   
    return expression_data, data['GeneSymbol'].values

# Performing 2-way ANOVA
def perform_Anova(expression_data):
    interaction_p_values = []
 
    # Here we define the design matrices
    intercept = np.ones(48)  # Intercept (Overall mean)
    smoking_status = np.array([0] * 12 + [1] * 12 + [0] * 12 + [1] * 12)  # 0 for Non-Smokers; 1 for Smokers
    gender = np.array([1] * 24 + [0] * 24)  # 1 for Males; 0 for Females
 
    # Additive Model (N) matrix
    N = np.vstack([intercept, smoking_status, gender]).T
 
    # Interaction Model (D) matrix
    interaction = smoking_status * gender  # Interaction term: 1 if both Smoking=1 and Gender=1
    D = np.vstack([intercept, smoking_status, gender, interaction]).T
 
    # Precomputing Projection Matrices(PM)
    P_N = N @ np.linalg.pinv(N)  # PM for Additive Model
    P_D = D @ np.linalg.pinv(D)  # PM for Interaction Model
    I = np.identity(48)
 
    # Calculating degrees of freedom using matrix ranks
    rank_N = np.linalg.matrix_rank(N)
    rank_D = np.linalg.matrix_rank(D)
    df_interaction = rank_D - rank_N  # DOF for interaction effect
    df_residual = 48 - rank_D  # DOF for residuals
 
    # Looping through each probe and calculating the interaction p-values.
    for i in range(expression_data.shape[0]):
        Y = expression_data[i, :].reshape(-1, 1)
 
        # Calculating sum of squares using projection matrices and trace properties
        SS_interaction_effect = np.trace((P_D - P_N) @ (Y @ Y.T))
        SS_residual = np.trace((I - P_D) @ (Y @ Y.T))
        # Calculating the F-statistic
        MS_interaction = SS_interaction_effect / df_interaction
        MS_residual = SS_residual / df_residual # These are Mean Squares.
        F_statistic = MS_interaction / MS_residual
        p_value = 1 - f.cdf(F_statistic, df_interaction, df_residual) # p value calculation
 
        interaction_p_values.append(p_value)
 
    return interaction_p_values
 
# Running ANOVA.
expression_data, gene_symbols = load_and_preprocess_data('Raw Data_GeneSpring.txt')  # REPLACE WITH YOUR APPROPRIATE FILEPATH!!!
p_values = perform_Anova(expression_data)
 
# Histogram plotting
plt.figure(figsize=(10, 6))
plt.hist(p_values, bins=100, edgecolor='black', color = 'orange')
plt.xlabel('Interaction p-values (Smoking Status x Gender)')
plt.ylabel('Frequency')
plt.title('Histogram of p-values for Gene Interaction')
plt.show()
