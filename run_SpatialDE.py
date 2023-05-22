import pandas as pd
import os
import sys

#rcParams['axes.spines.right'] = False
#rcParams['axes.spines.top'] = False

import NaiveDE
import SpatialDE


for data_n in [20] + list(range(61, 81)) + [24]:
    try:
        print('Data', data_n)

        folder = '../SPIRAL_webtool/static/analysis/data' + str(data_n)

        out_folder = 'comparison_to_SpatialDE/data' + str(data_n) + '/'
        if not os.path.isdir('./comparison_to_SpatialDE'):
            os.mkdir('./comparison_to_SpatialDE')
        if not os.path.isdir(out_folder):
            os.mkdir(out_folder)

        # load raw file
        counts_file = os.path.join(folder, 'counts.csv')
        counts = pd.read_csv(counts_file, index_col=0)
        counts = counts.loc[counts.sum(axis=1) != 0]
        counts = counts.groupby(level=0).sum()

        total_counts = counts.sum()

        counts = counts[counts.sum(1) >= 3].T  # Filter practically unobserved genes

        pos_file = os.path.join(folder, 'spatial_coors.csv')
        if data_n==5:
            pos = pd.read_csv(pos_file, index_col=0, sep='\t', header=None)
        else:
            pos = pd.read_csv(pos_file, index_col=0)
        pos.columns = ['x', 'y']
        pos.index.name = None
        pos = pos.loc[counts.index, :]
        pos.loc[:, 'total_counts'] = total_counts

        #figsize(6, 4)
        #plt.scatter(pos['x'], pos['y'], c='k');
        #plt.axis('equal');

        norm_expr = NaiveDE.stabilize(counts.T).T
        resid_expr = NaiveDE.regress_out(pos, norm_expr.T, 'np.log(total_counts)').T

        X = pos[['x', 'y']]
        results = SpatialDE.run(X, resid_expr)
        sign_results = results.query('qval < 0.05')

        for C in [3, 5]:
            histology_results, patterns = SpatialDE.aeh.spatial_patterns(X, resid_expr, sign_results, C=C, l=1.8, verbosity=1)

            pattern_list = list(set(histology_results['pattern']))
            method_cluster_table = pd.DataFrame(columns=['#genes', 'genes'], index=pattern_list)
            for m in pattern_list:
                genes = list(histology_results[histology_results['pattern'] == m].g)
                method_cluster_table.loc[m, "#genes"] = len(genes)
                method_cluster_table.loc[m, "genes"] = str(genes)
            method_cluster_table.to_csv(os.path.join(out_folder, 'cluster_table_' + str(C) + '.csv'))

        del counts

    except Exception as e:
        out_folder = 'comparison_to_SpatialDE/data' + str(data_n) + '/'
        if not os.path.isdir('./comparison_to_SpatialDE'):
            os.mkdir('./comparison_to_SpatialDE')
        if not os.path.isdir(out_folder):
            os.mkdir(out_folder)
        with open(os.path.join(out_folder, "error.txt"), "w") as text_file:
            text_file.write(str(e))

