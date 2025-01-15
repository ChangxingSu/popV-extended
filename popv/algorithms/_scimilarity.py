from __future__ import annotations

import logging
import scanpy as sc
import numpy as np
from scimilarity import CellAnnotation # import CellAnnotation from scimilarity
from scimilarity.utils import lognorm_counts, align_dataset # import lognorm_counts and align_dataset from scimilarity.utils

class SCIMILARITY_POPV:
    def __init__(
        self,
        batch_key: str | None = "_batch_annotation",
        labels_key: str | None = "_labels_annotation", 
        model_path: str | None = None,
        result_key: str | None = "popv_scimilarity_prediction",
        embedding_key: str | None = "X_scimilarity_umap_popv",
        target_celltypes: list | None = None,
    ) -> None:
        """
        Class to compute cell type predictions using SCimilarity.

        Parameters
        ----------
        batch_key
            Key in obs field of adata for batch information.
        labels_key
            Key in obs field of adata for cell-type information.
        model_path
            Path to SCimilarity pretrained model.
        result_key
            Key in obs in which celltype annotation results are stored.
        embedding_key
            Key in obsm in which UMAP embedding is stored.
        target_celltypes
            List of target cell types to constrain classification.
        method_dict
            Additional parameters for SCimilarity.
        embedding_dict
            Dictionary for UMAP parameters.
        """
        self.batch_key = batch_key
        self.labels_key = labels_key
        self.model_path = model_path
        self.result_key = result_key
        self.embedding_key = embedding_key
        self.target_celltypes = target_celltypes

    def compute_integration(self, adata):
        """
        No integration needed, just compute embeddings using pretrained model
        """
        logging.info("Computing SCimilarity embeddings using pretrained model")
        
        # Initialize SCimilarity model
        self.model = CellAnnotation(model_path=self.model_path)
        
        # Preprocess data
        adata_aligned = align_dataset(adata, self.model.gene_order)
        adata_aligned = lognorm_counts(adata_aligned)
        
        # Compute embeddings using pretrained model
        adata.obsm["X_scimilarity"] = self.model.get_embeddings(adata_aligned.X)

        if self.target_celltypes is not None:
            self.model.safelist_celltypes(self.target_celltypes)

    def predict(self, adata):
        """
        Make predictions using pretrained model
        """
        logging.info(
            f'Computing SCimilarity predictions using pretrained model. Storing in adata.obs["{self.result_key}"]'
        )

        # Get predictions using pretrained model
        predictions, nn_idxs, nn_dists, nn_stats = self.model.get_predictions_knn(
            adata.obsm["X_scimilarity"]
        )
        
        # Store results
        adata.obs[self.result_key] = predictions.values
        
        # Store QC metrics
        adata.obs[self.result_key + "_min_dist"] = nn_dists
        
        if adata.uns["_return_probabilities"]:
            # Store prediction probabilities from nn_stats
            adata.obs[self.result_key + "_probabilities"] = nn_stats["confidence"].values

    def compute_embedding(self, adata):
        """
        Compute UMAP from SCimilarity embeddings
        """
        if adata.uns["_compute_embedding"]:
            logging.info(
                f'Computing UMAP from SCimilarity embeddings. Storing in adata.obsm["{self.embedding_key}"]'
            )
            sc.pp.neighbors(adata, use_rep="X_scimilarity")
            adata.obsm[self.embedding_key] = sc.tl.umap(
                adata, copy=True
            ).obsm["X_umap"] 