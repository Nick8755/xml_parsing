import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Creating a DataFrame to summarize the relationship between cancer cell lines and bazedoxifene concentrations
data = {
    "Cancer Cell Line": [
        "DAOY (Medulloblastoma)",
        "UW288 (Medulloblastoma)",
        "MDA-MB-231 (Triple-Negative Breast Cancer)",
        "PANC-1 (Pancreatic Cancer)"
    ],
    "Effective Concentration (µM)": [
        "20 µM (Inhibition of IL-6 signaling)",
        "5.65 µM (IC50 for proliferation)",
        "10-20 µM (Inhibition of IL-6 signaling)",
        "20 µM (Inhibition of IL-6 signaling)"
    ]
}

cancer_concentration_table = pd.DataFrame(data)
cancer_concentration_table

#print(cancer_concentration_table)



# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

# Draw the components of the IL-6 signaling pathway
# IL-6
ax.annotate('IL-6', xy=(0.1, 0.8), fontsize=12, ha='center', color='blue')
ax.add_patch(mpatches.Circle((0.1, 0.8), 0.05, color='lightblue', ec='blue'))

# IL-6R
ax.annotate('IL-6R', xy=(0.1, 0.6), fontsize=12, ha='center', color='blue')
ax.add_patch(mpatches.Circle((0.1, 0.6), 0.05, color='lightblue', ec='blue'))

# GP130
ax.annotate('GP130', xy=(0.1, 0.4), fontsize=12, ha='center', color='blue')
ax.add_patch(mpatches.Circle((0.1, 0.4), 0.05, color='lightblue', ec='blue'))

# STAT3
ax.annotate('STAT3', xy=(0.1, 0.2), fontsize=12, ha='center', color='blue')
ax.add_patch(mpatches.Circle((0.1, 0.2), 0.05, color='lightblue', ec='blue'))

# Bazedoxifene
ax.annotate('Bazedoxifene', xy=(0.6, 0.6), fontsize=12, ha='center', color='red')
ax.add_patch(mpatches.Rectangle((0.55, 0.55), 0.1, 0.1, color='salmon', ec='red'))

# Arrows indicating action
ax.annotate('', xy=(0.1, 0.75), xytext=(0.1, 0.65), arrowprops=dict(arrowstyle='<-', color='black'))
ax.annotate('', xy=(0.1, 0.55), xytext=(0.1, 0.45), arrowprops=dict(arrowstyle='->', color='black'))
ax.annotate('', xy=(0.1, 0.35), xytext=(0.1, 0.25), arrowprops=dict(arrowstyle='->', color='black'))

# Indicate inhibition
ax.annotate('Inhibition of IL-6/GP130\nsignaling pathway', xy=(0.6, 0.4), fontsize=12, ha='center', color='red')
ax.annotate('', xy=(0.6, 0.55), xytext=(0.6, 0.45), arrowprops=dict(arrowstyle='->', color='red'))

# Set limits and remove axes
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

# Title
plt.title('Action of Bazedoxifene on IL-6 Signaling Pathway', fontsize=14)

# Show the plot
plt.show()