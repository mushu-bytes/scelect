import sys
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import altair as alt
from anndata import AnnData
from typing import List
from scprocessing.SelectPipeline import SelectPipeline
from utils import splitAD, read_single_cell_data


if __name__ == "__main__":
    # arguments would be input file, output file, key
    input_path, output_path, dataset_key = sys.argv[1], sys.argv[2], sys.argv[3]
    immune = read_single_cell_data(f"{input_path}/human_brca_immune.h5ad")
    immune.var_names_make_unique()
    immune.obs_names_make_unique()
    del immune.obsm["X_diffmap"] # if available throws off script
    immune_data = splitAD(immune, dataset_key)[0:2]
    select_immune = SelectPipeline(normalization=["seurat", "zheng17"],\
                                    integration=["harmony", "scanorama", "merge"],\
                                    metrics=["jaccard", "ari", "nmi"],
                                    resolution_range=[0.3],
                                    key=dataset_key
                                )
    immune_res, report_immune, pipeline_immune = select_immune.search(immune_data, key_metric="ari")
    immune_res.write(f"{output_path}/scanpy_integration.h5ad")

    # TODO: Move visualization to separate file
    # post processing for visualization, probably move this into a different file
    final_df = []
    for key, cluster in select_immune.clusters.items():
        init_df = pd.DataFrame(cluster.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
        init_df[dataset_key] = cluster.obs.tissue_condition.tolist()
        init_df["Normalization"] = key[0]
        init_df["Integration"] = key[1]
        final_df.append(init_df)
    final_df = pd.concat(final_df)

    # this cell just filters the final_df
    all_cells = []
    for condition in final_df[dataset_key].unique():
        all_cells.append(final_df[final_df[dataset_key] == condition].sample(3000))
    final_df = pd.concat(all_cells)

    alt.data_transformers.enable("vegafusion")
    high_contrast_colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231']

    umap_chart = alt.Chart(final_df).mark_point(size=1.5).encode(
        x="UMAP1",
        y="UMAP2",
        color=alt.Color(dataset_key).scale(range=high_contrast_colors)
    )

    facet_grid = umap_chart.facet(
        row="Normalization",
        column="Integration",
        title='UMAP plots for Different Combinations of Normalization and Integration Methods'
    )
    facet_grid = facet_grid.configure_axis(
        grid=False,
        labelFontSize=18,
        titleFontSize=18,
        tickSize=2
    )

    facet_grid = facet_grid.configure_header(
        labelFontSize=18,
        titleFontSize=18
    )

    facet_grid = facet_grid.configure_headerRow(
        labelFontSize=18,
        titleFontSize=18
    )

    facet_grid = facet_grid.configure_headerColumn(
        labelFontSize=18,
        titleFontSize=18
    )

    facet_grid = facet_grid.configure_legend(
        titleFontSize=18,
        labelFontSize=16
    )

    facet_grid.save(f"{output_path}/scanpy_report.png", scale_factor=5.0)
    