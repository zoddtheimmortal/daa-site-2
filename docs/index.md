# Efficient Algorithms for Densest Subgraph Discovery

## DSD Problem

Densest subgraph discovery (DSD) is a fundamental problem in graph mining that has been studied for decades. It finds widespread applications across various domains including network science, biological analysis, and graph databases. The core objective of DSD is to identify a subgraph D within a given graph G that exhibits the highest density (typically measured as the ratio of edges to vertices in D).

Due to the inherent complexity of solving DSD, we introduce a novel solution paradigm in this website. Our key insight reveals that the densest subgraph can be accurately identified through a k-core (a specific type of dense subgraph of G), accompanied by theoretical guarantees. Building on this intuition, we develop both efficient exact and approximation solutions for addressing the DSD problem.

A notable strength of our proposed solutions is their versatility in finding densest subgraphs across a broad spectrum of graph density definitions, encompassing both clique-based and general pattern-based density metrics. We have conducted comprehensive experimental evaluations using both real-world and synthetic datasets. The results demonstrate that our algorithms achieve performance improvements of up to four orders of magnitude compared to existing approaches.

## Algorithms

-   [Algorithm 1 - Exact](algo1.md)
-   [Algorithm 4 - CoreExact](algo4.md)
