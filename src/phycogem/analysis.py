import matplotlib.pyplot as plt
import seaborn as sns
from pandas import DataFrame


def plot_flux_distribution(
    sample: DataFrame,
    reaction_ids: list,
    histogram: bool = False,
    fva: DataFrame = None,
    figsize_per_plot: tuple = (10, 6),
):
    """Plot the flux ranges and sample for a list of reactions.

    Args:
        sample (DataFrame): DataFrame with sampled flux distributions.
        reaction_ids (list): List of reaction IDs to plot.
        histogram (bool, optional): Whether to overlay a histogram. Defaults to False.
        fva (DataFrame, optional): DataFrame with flux variability analysis results. Defaults to None.
        figsize_per_plot (tuple, optional): Size of each individual plot. Defaults to (10, 6).
    """

    n_reactions = len(reaction_ids)
    n_rows = (n_reactions + 1) // 2

    fig, axes = plt.subplots(
        n_rows, 2, figsize=(figsize_per_plot[0] * 2, figsize_per_plot[1] * n_rows)
    )
    if n_reactions == 1:
        axes = [[axes]]
    elif n_reactions % 2 == 1:
        fig.delaxes(axes[-1, -1])

    for i, reaction_id in enumerate(reaction_ids):
        ax = axes[i // 2, i % 2]

        if histogram:
            sns.histplot(
                sample[reaction_id],
                kde=True,
                bins=30,
                color="skyblue",
                label="Histogram",
                ax=ax,
            )
        sns.kdeplot(
            sample[reaction_id],
            color="blue",
            label="Density Estimation",
            fill=True,
            ax=ax,
        )

        if fva is not None:
            min_flux = fva.loc[reaction_id, "minimum"]
            max_flux = fva.loc[reaction_id, "maximum"]
            ax.axvline(
                min_flux, color="r", linestyle="--", label=f"Min Flux: {min_flux:.2f}"
            )
            ax.axvline(
                max_flux, color="b", linestyle="--", label=f"Max Flux: {max_flux:.2f}"
            )

        ax.set_title(f"Flux Distribution for {reaction_id}")
        ax.set_xlabel("mmol/gDW/h")
        ax.set_ylabel("Density")
        ax.legend()

    plt.tight_layout()
    plt.show()
