import matplotlib.pyplot as plt
import pandas as pd
from numpy import arange

data = pd.read_csv('solv_data.csv')
print(data.columns)

X_positions = arange(0, len(list(data.structure)))

fig, ax = plt.subplots()
# ax.scatter(X_positions, data.gas-data.gas, c='k', marker='o',
#            edgecolor='k', s=80, label='APFD - gas phase')
# ax.scatter(X_positions, data.gas_gfn-data.gas_gfn, c='r', marker='X',
#            edgecolor='k', s=80, label='xTB - gas phase')
ax.scatter(X_positions, data.dcm_pcm-data.gas, c='k', marker='o',
           edgecolor='k', s=80, label='PCM - DCM')
ax.scatter(X_positions, data.dcm_cos-data.gas, c='r', marker='D',
           edgecolor='k', s=80, label='COSMO - DCM')
ax.scatter(X_positions, data.dcm_gfn-data.gas_gfn, c='b', marker='X',
           edgecolor='k', s=80, label='xTB - DCM')
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlabel('structure', fontsize=16)
ax.set_ylabel('free energy of formation \n (solvated - gas) [kJ/mol]',
              fontsize=16)
ax.axhline(y=0, c='k', alpha=0.2)
# ax.set_xlim(0, 8)
ax.set_ylim(-250, 50)
ax.legend()
ax.set_xticks(X_positions)
ax.set_xticklabels(X_positions)
fig.tight_layout()
fig.savefig('solvation_comparison_dcm.pdf', dpi=720, bbox_inches='tight')

fig, ax = plt.subplots()
# ax.scatter(X_positions, data.gas-data.gas, c='k', marker='o',
#            edgecolor='k', s=80, label='APFD - gas phase')
# ax.scatter(X_positions, data.gas_gfn-data.gas_gfn, c='r', marker='X',
#            edgecolor='k', s=80, label='xTB - gas phase')
ax.scatter(X_positions, data.chcl3_pcm-data.gas, c='k', marker='o',
           edgecolor='k', s=80, label='PCM - CHCl$_3$')
ax.scatter(X_positions, data.chcl3_cos-data.gas, c='r', marker='D',
           edgecolor='k', s=80, label='COSMO - CHCl$_3$')
ax.scatter(X_positions, data.chcl3_gfn-data.gas_gfn, c='b', marker='X',
           edgecolor='k', s=80, label='xTB - CHCl$_3$')
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlabel('structure', fontsize=16)
ax.set_ylabel('free energy of formation \n (solvated - gas) [kJ/mol]',
              fontsize=16)
ax.axhline(y=0, c='k', alpha=0.2)
# ax.set_xlim(0, 8)
ax.set_ylim(-250, 50)
ax.legend()
ax.set_xticks(X_positions)
ax.set_xticklabels(X_positions)
fig.tight_layout()
fig.savefig('solvation_comparison_chcl3.pdf', dpi=720, bbox_inches='tight')

fig, ax = plt.subplots()
# ax.scatter(X_positions, data.gas-data.gas, c='k', marker='o',
#            edgecolor='k', s=80, label='APFD - gas phase')
# ax.scatter(X_positions, data.gas_gfn-data.gas_gfn, c='r', marker='X',
#            edgecolor='k', s=80, label='xTB - gas phase')
ax.scatter(X_positions, data.thf_pcm-data.gas, c='k', marker='o',
           edgecolor='k', s=80, label='PCM - THF')
ax.scatter(X_positions, data.thf_cos-data.gas, c='r', marker='D',
           edgecolor='k', s=80, label='COSMO - THF')
ax.scatter(X_positions, data.thf_gfn-data.gas_gfn, c='b', marker='X',
           edgecolor='k', s=80, label='xTB - THF')
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlabel('structure', fontsize=16)
ax.set_ylabel('free energy of formation \n (solvated - gas) [kJ/mol]',
              fontsize=16)
ax.axhline(y=0, c='k', alpha=0.2)
# ax.set_xlim(0, 8)
ax.set_ylim(-250, 50)
ax.legend()
ax.set_xticks(X_positions)
ax.set_xticklabels(X_positions)
fig.tight_layout()
fig.savefig('solvation_comparison_thf.pdf', dpi=720, bbox_inches='tight')
