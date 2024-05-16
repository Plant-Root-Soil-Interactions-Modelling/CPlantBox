from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from statannotations.Annotator import Annotator
from scikit_posthocs import posthoc_tukey
import seaborn as sns
from scipy.stats import linregress
import os

dir = os.path.dirname(os.path.abspath(__file__))

'''FIGURE 3'''
# 'load data'
idx = ['P0','P1','P2','P3']
a_axial = pd.read_csv(dir + '/a_axial.csv')
r_crown = pd.read_csv(dir + '/r_crown.csv')
r_leaves = pd.read_csv(dir + '/r_leaves.csv')
rootshoot = pd.read_csv(dir + '/root_shoot_ratio.csv')

'''PANEL A'''
plt.subplots(figsize=(6,4))
tukey_df = posthoc_tukey(a_axial, val_col="value", group_col="variable")
remove = np.tril(np.ones(tukey_df.shape), k=0).astype("bool")
tukey_df[remove] = np.nan
molten_df = tukey_df.melt(ignore_index=False).reset_index().dropna()
ax = sns.boxplot(data=a_axial, x="variable", y="value", order=idx)
pairs = [(i[1]["index"], i[1]["variable"]) for i in molten_df.iterrows()]
p_values = [i[1]["value"] for i in molten_df.iterrows()]
annotator = Annotator(
    ax, pairs, data=a_axial, x="variable", y="value", order=idx)
annotator.configure(text_format="star", loc="inside")
annotator.set_pvalues_and_annotate(p_values)
plt.xlabel('P level')
plt.ylabel('axial root diameter (cm)')
plt.tight_layout()
plt.show()

'''PANEL B'''
plt.subplots(figsize=(6,4))
tukey_df = posthoc_tukey(r_crown, val_col="value", group_col="variable")
remove = np.tril(np.ones(tukey_df.shape), k=0).astype("bool")
tukey_df[remove] = np.nan
molten_df = tukey_df.melt(ignore_index=False).reset_index().dropna()
ax = sns.boxplot(data=r_crown, x="variable", y="value", order=idx)
pairs = [(i[1]["index"], i[1]["variable"]) for i in molten_df.iterrows()]
p_values = [i[1]["value"] for i in molten_df.iterrows()]
annotator = Annotator(
    ax, pairs, data=r_crown, x="variable", y="value", order=idx)
annotator.configure(text_format="star", loc="inside")
annotator.set_pvalues_and_annotate(p_values)
plt.xlabel('P level')
plt.ylabel('Crown r (cm d-1)')
plt.tight_layout()
plt.show()

'''PANEL C'''
plt.subplots(figsize=(6,4))
tukey_df = posthoc_tukey(r_leaves, val_col="value", group_col="variable")
remove = np.tril(np.ones(tukey_df.shape), k=0).astype("bool")
tukey_df[remove] = np.nan
molten_df = tukey_df.melt(ignore_index=False).reset_index().dropna()
ax = sns.boxplot(data=r_leaves, x="variable", y="value", order=idx)
pairs = [(i[1]["index"], i[1]["variable"]) for i in molten_df.iterrows()]
p_values = [i[1]["value"] for i in molten_df.iterrows()]
annotator = Annotator(
    ax, pairs, data=r_leaves, x="variable", y="value", order=idx)
annotator.configure(text_format="star", loc="inside")
annotator.set_pvalues_and_annotate(p_values)
plt.xlabel('P level')
plt.ylabel('leaf elongation rate (cm day$^{-1}$)')
plt.tight_layout()
plt.show()

'''PANEL D'''
rootshoot['root_shoot'] = rootshoot.root_dm / rootshoot.shoot_dm
grouped = rootshoot.groupby('level')
def standard_error(x):
    return np.std(x, ddof=1) / np.sqrt(x.count())
result = grouped.agg({'shoot_dm': ['mean', standard_error], 'root_dm': ['mean', standard_error], 'root_shoot': ['mean', standard_error]})
result.columns = ['shoot_dm_mean', 'shoot_dm_sem', 'root_dm_mean', 'root_dm_sem', 'root_shoot_mean', 'root_shoot_sem']
result.reset_index(inplace=True)
x_values = [1.8, 3.2, 4.4, 7.7]
slope, intercept, r_value, p_value, std_err = linregress(x_values, result['root_shoot_mean'])
y_fit = intercept + slope * np.array(x_values)

plt.figure(figsize=(10, 6))
ax = plt.subplot()
plt.errorbar(
    x_values, result['root_shoot_mean'], 
    yerr=result['root_shoot_sem'], 
    fmt='o', label='root:shoot ratio')
for x, level in zip(x_values, ['P0', 'P1', 'P2', 'P3']):
    plt.text(x, result['root_shoot_mean'].loc[result['level'] == level].values[0] + 0.1, level, ha='center')
plt.plot(x_values, y_fit, linestyle='--', color='red', label='Linear Regression Line')
plt.text(
    2.5, 1.3, f"R²: {r_value**2:.2f}", 
    bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
ax.set_xlabel('Pcal mg 100g$^{-1}$')
ax.set_ylabel('[-]')
plt.xlim(1.5, 8)
plt.ylim(0.6, 2)
plt.legend()
plt.tight_layout()
plt.show()


'''FIGURE 5'''
'''PLOT PANEL A'''
data = pd.read_excel('C:/Users/mobil/Desktop/Experiment_Rhizotrons/response_curves/response_data.xlsx')
adjusted_data = data[(data['Pcal'] >= 1.8) & (data['Pcal'] <= 7.3)]
X = adjusted_data[['Pcal']]
y = adjusted_data['a_axial']
X_with_const = sm.add_constant(X)
model = sm.OLS(y, X_with_const).fit()
print(model.summary())
p_value = model.pvalues[1]
print(f"P-value for 'Pcal': {p_value}")
pcal_values = np.linspace(1.8, 7.3, 500)
X_plot = sm.add_constant(pcal_values)  # Remember to add constant for prediction
predictions = model.predict(X_plot)
r_squared = model.rsquared

plt.figure(figsize=(10, 6))
plt.scatter(X, y, color='gray', alpha=0.5, label='Data')
plt.plot(pcal_values, predictions, color='black', linewidth=2)
plt.text(x=0.5, y=0.1, s=f'R² = {r_squared:.4f}\nP-value = {p_value:.4g}', fontsize=12, transform=plt.gca().transAxes)
plt.xlabel('P$_{Cal}$ (g 100g$^{-1}$)')
plt.ylabel('axial root radii (cm)')
plt.xlim(1.6, 7.5)
plt.ylim(0.02, 0.12)
plt.show()

'''PLOT PANEL B'''
X = np.array(data['P_DM']).reshape(-1, 1)
y = np.array(data['r_crown']).reshape(-1, 1)
model_origin = LinearRegression(fit_intercept=False)
model_origin.fit(X, y)
y_pred_origin = model_origin.predict(X)
r2_origin = r2_score(y, y_pred_origin)
coef_origin = model_origin.coef_[0]
intercept_origin = 0  # Explicitly set to 0
X_plot = np.array([[0], [X.max()]])  # From 0 to max of X
y_pred_plot = model_origin.predict(X_plot)  # Predict using the model
print(r2_origin, coef_origin, intercept_origin)

plt.figure(figsize=(10, 6))
plt.scatter(X, y, color='grey')  # Data points
plt.plot(X_plot, y_pred_plot, color='black', linewidth=2)  # Plotting the line
plt.xlabel('P$_{Cal}$/DM')
plt.ylabel('crown root elongation (cm d^-1)')
plt.xlim(0, 0.0011)
plt.ylim(0, 7.8)
plt.ticklabel_format(axis= 'both',style = 'scientific',scilimits=(-2,5))
plt.show()

'''PANEL C'''
new_model_data = data[['Pcal', 'r_leaf']]
X_Pcal_with_const = sm.add_constant(new_model_data[['Pcal']])
model = sm.OLS(new_model_data['r_leaf'], X_Pcal_with_const).fit()
r_squared = model.rsquared
p_values = model.pvalues
print(f"R-squared: {r_squared}")
print(f"P-values:\n{p_values}")
Pcal_observed_range = np.linspace(new_model_data['Pcal'].min(), new_model_data['Pcal'].max(), 100)
r_leaf_pred = model.predict(sm.add_constant(Pcal_observed_range))

plt.figure(figsize=(10, 6))
plt.scatter(new_model_data['Pcal'], new_model_data['r_leaf'], color='grey')
plt.plot(Pcal_observed_range, r_leaf_pred, color='black')
plt.xlabel('P$_{Cal}$ (g 100g$^{-1}$)')
plt.ylabel('leaf elongation rate (cm d^-1)$')
plt.xlim(1.6, 7.5)
plt.ylim(2.5, 15.5)
plt.show()
