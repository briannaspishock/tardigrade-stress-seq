import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Set figure and canvas
fig, ax = plt.subplots(figsize=(11, 6.5), dpi=300)
ax.set_facecolor('#0d1117')
fig.patch.set_facecolor('#0d1117')

# Draw Background Compartments
ax.axhspan(5.5, 9.5, facecolor='#161b22', alpha=0.9)
ax.axhspan(0.5, 4.5, facecolor='#0f1c24', alpha=0.9)

# Draw Plasma Membrane
membrane = patches.Rectangle((0.5, 4.5), 9.0, 1.0, linewidth=1.5, edgecolor='#58a6ff', facecolor='#1f2937', zorder=2)
ax.add_patch(membrane)
ax.text(5.0, 5.0, "PLASMA MEMBRANE (LIPID BILAYER)", color='#58a6ff', fontsize=11, fontweight='bold', ha='center', va='center', zorder=3, family='sans-serif')

# === 1. Extracellular / SAHS Compartment (Top) ===
ax.text(0.8, 9.0, "EXTRACELLULAR SPACE / HEMOCOEL", color='#8b949e', fontsize=9.5, fontweight='bold', family='sans-serif', va='top')
ax.text(0.8, 8.4, "SAHS: Secretory Abundant Heat Soluble", color='#3fb950', fontsize=12, fontweight='bold', family='sans-serif', va='top')
sahs_bullets = "• Contains N-terminal secretion signal\n• Secreted externally to shield membrane leaflets\n• Prevents desiccation-induced membrane fusion"
ax.text(0.8, 7.7, sahs_bullets, color='#c9d1d9', fontsize=9, linespacing=1.5, va='top')

# Draw SAHS Shielding Particles
np.random.seed(42)
sahs_x = np.random.uniform(5.5, 9.0, 18)
sahs_y = np.random.uniform(5.7, 7.4, 18)
ax.scatter(sahs_x, sahs_y, c='#3fb950', s=90, edgecolors='#ffffff', linewidth=0.8, alpha=0.9, zorder=4, label='SAHS (Secreted Peptides)')

# Secretion Arrow
ax.annotate('', xy=(7.2, 5.7), xytext=(7.2, 4.2),
            arrowprops=dict(facecolor='#3fb950', edgecolor='#3fb950', width=2, headwidth=6), zorder=5)
ax.text(7.2, 3.9, "Secretion Outward", color='#3fb950', fontsize=8.5, fontweight='bold', ha='center', zorder=5)

# === 2. Cytosolic / CAHS Compartment (Bottom) ===
ax.text(0.8, 4.0, "CYTOPLASM (INTRACELLULAR)", color='#8b949e', fontsize=9.5, fontweight='bold', family='sans-serif', va='top')
ax.text(0.8, 3.4, "CAHS: Cytoplasmic Abundant Heat Soluble", color='#f778ba', fontsize=12, fontweight='bold', family='sans-serif', va='top')
cahs_bullets = "• Lacks secretion signal (cytosolic retention)\n• Intrinsically disordered coils in hydration\n• Reversible self-assembly into protective Bioglass matrix\n• Vitrifies to entrap and protect native enzymes"
ax.text(0.8, 2.7, cahs_bullets, color='#c9d1d9', fontsize=9, linespacing=1.5, va='top')

# Draw CAHS Network
cahs_x = np.random.uniform(5.5, 9.0, 16)
cahs_y = np.random.uniform(1.2, 3.5, 16)
for i in range(len(cahs_x)-1):
    for j in range(i+1, len(cahs_x)):
        if np.hypot(cahs_x[i]-cahs_x[j], cahs_y[i]-cahs_y[j]) < 1.6:
            ax.plot([cahs_x[i], cahs_x[j]], [cahs_y[i], cahs_y[j]], color='#f778ba', alpha=0.35, lw=1.5, zorder=3)
ax.scatter(cahs_x, cahs_y, c='#f778ba', s=95, edgecolors='#ffffff', linewidth=0.8, alpha=0.9, zorder=4, label='CAHS (Cytosolic Vitrification Network)')

# Entrapped Native Enzyme
ax.scatter([7.4], [2.4], c='#e3b341', s=220, marker='s', edgecolors='#ffffff', linewidth=1.2, zorder=5, label='Native Enzymes (Protected from Unfolding)')
ax.text(7.4, 1.85, "Entrapped Enzyme", color='#e3b341', fontsize=8, ha='center', fontweight='bold', zorder=5)

# Coordinate Boundaries & Title
ax.set_xlim(0.5, 9.5)
ax.set_ylim(0.5, 9.5)
ax.axis('off')

fig.suptitle("Dual-Compartment Anhydrobiosis Defense Architecture", fontsize=14, fontweight='bold', color='#f0f6fc', y=0.96)

# Legend
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.04), ncol=3, frameon=True, facecolor='#161b22', edgecolor='#30363d', fontsize=8.5)
for text in leg.get_texts():
    text.set_color('#c9d1d9')

plt.tight_layout()
plt.savefig("figures/cahs_vs_sahs_mechanism.png", dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
print("[+] Generated clean mechanism figure: figures/cahs_vs_sahs_mechanism.png")
