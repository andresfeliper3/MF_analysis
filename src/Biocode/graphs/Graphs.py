import os
from itertools import cycle

import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
import pandas as pd
import seaborn as sns

from utils.logger import logger


class Graphs:

    @staticmethod
    def _savefig(title, name):
        directory = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), "out/graphs")
        os.makedirs(directory, exist_ok=True)

        actual_path = f"{directory}/{name}"
        plt.tight_layout()
        if not os.path.exists(actual_path):
            os.makedirs(actual_path)
        plt.ioff()
        plt.savefig(f'{actual_path}/{title}.png')
        plt.close()

    @staticmethod
    def graph_one(x_array, y_array, x_label, y_label, title, name, save=True):
        plt.figure(figsize=(10, 6))
        plt.plot(x_array, y_array, marker='o')
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.grid()
        if save:
            Graphs._savefig(title, name)
        plt.show()

    @staticmethod
    def graph_many(results_array, X, Y, x_label, y_label, title, name, markers_array=None, linestyles_array=None,
                   colors_array=None, labels_array=None, markersize=6, save=True):
        plt.figure(figsize=(10, 6))

        if markers_array is None:
            markers = ['o', 's', '^', 'v', '>', '<', 'p', 'D', 'h']
            marker_cycle = cycle(markers)
        else:
            marker_cycle = cycle(markers_array)

        for index, result in enumerate(results_array):
            marker = next(marker_cycle)
            linestyle = linestyles_array[index] if linestyles_array else '-'
            color = colors_array[index] if colors_array else None
            label = labels_array[index] if labels_array else None

            plt.plot(result[X], result[Y], marker=marker, linestyle=linestyle, color=color, label=label,
                     markersize=markersize)

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.grid()
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1),  ncol=2)
        if save:
            Graphs._savefig(title, name)
        plt.show()

    @staticmethod
    def graph_many_grouped(results_array, X, Y, x_label, y_label, title, name, regions_number=None, markers_array=None,
                           linestyles_array=None,
                           colors_array=None, labels_array=None, markersize=6, color_by='region', save=True):
        if not (regions_number >= 0):
            raise Exception("Not a valid regions_number entered in the graph_many_grouped method of Graphs")

        plt.figure(figsize=(10, 6))

        markers = ['o', 's', '^', 'v', '>', '<', 'p', 'D', 'h']
        colors = ['b', 'r', 'g', 'c', 'y', 'm', 'k', 'w']

        markers_cycle = cycle(markers[:regions_number]) if markers_array is None else cycle(markers_array)
        colors_cycle = cycle(colors[:regions_number]) if colors_array is None else cycle(colors_array)

        for index, result in enumerate(results_array):
            if color_by == 'region':
                if index % regions_number == 0:
                    marker = markers_array[index] if markers_array else next(markers_cycle)
                color = colors_array[index] if colors_array else next(colors_cycle)
            elif color_by == 'chromosome':
                if index % regions_number == 0:
                    color = colors_array[index] if colors_array else next(colors_cycle)
                marker = markers_array[index] if markers_array else next(markers_cycle)

            linestyle = linestyles_array[index] if linestyles_array else '-'
            label = labels_array[index] if labels_array else None
            plt.plot(result[X], result[Y], marker=marker, linestyle=linestyle, color=color, label=label,
                     markersize=markersize)

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.grid()
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1),  ncol=2)
        if save:
            Graphs._savefig(title, name)
        plt.show()

    @staticmethod
    def graph_bars(x_array, y_array, title, name, y_label=None, bar_labels=None, bar_colors=None, legend=None, rotation=45,
                   y_range: list[int] = None,
                   top_labels=False, save=True):
        fig, ax = plt.subplots()
        ax.bar(x_array, y_array, label=bar_labels, color=bar_colors)
        ax.set_ylabel(y_label)

        if y_range:
            ax.set_ylim(y_range[0], y_range[1])
        else:
            # Calculate the dynamic ylim based on the values in y_array
            min_y = min(y_array)
            max_y = max(y_array)
            buffer = 0.1
            ax.set_ylim(min_y - buffer, max_y + buffer)

        ax.set_title(title)
        # Calculate the x-positions and y-values for the labels
        if top_labels:
            for x, y in zip(x_array, y_array):
                ax.text(x, y, f"{y:.2f}", ha="center", va="bottom")

        ax.set_xticks(x_array)
        ax.set_xticklabels(x_array, rotation=rotation)
        ax.legend(title=legend)
        if save:
            Graphs._savefig(title, name)
        plt.show()

    @staticmethod
    def graph_bars_grouped(x_array, y_array, title, name, regions_number=3, y_label=None, x_labels=None,
                           legend_labels=None,
                           regions_colors=None, rotation=45, y_range: list[int] = None, top_labels=False, save=True):
        X_SIZE = 20
        fig, ax = plt.subplots(figsize=(X_SIZE, 6))  # Adjust figure size for better visibility
        bar_width = 0.1  # Width of each bar
        num_chromosomes = len(x_array) // regions_number

        if x_labels is None:
            x_labels = [f'Chromosome {i + 1}' for i in range(num_chromosomes)]

        # Generate a colormap based on the regions_number
        cmap = cm.get_cmap('viridis', regions_number)
        colors = [cmap(i / regions_number) for i in range(regions_number)]

        # Calculate the x-positions for the bars
        x_positions = [i * (regions_number + 2) + np.arange(regions_number) for i in range(num_chromosomes)]

        for i in range(num_chromosomes):
            for j in range(regions_number):
                color = colors[j]
                x_position = x_positions[i][j]
                y_value = y_array[i * regions_number + j]
                label = x_labels[i]

                ax.bar(x_position, y_value, label=label, color=color)
                # labels on top of each bar
                if top_labels:
                    ax.text(x_position, y_value, f"{y_value:.3f}", ha="center", va="bottom")

        ax.set_ylabel(y_label)

        if y_range:
            ax.set_ylim(y_range[0], y_range[1])
        else:
            # Calculate the dynamic ylim based on the values in y_array
            min_y = min(y_array)
            max_y = max(y_array)
            buffer = 0.1
            ax.set_ylim(min_y - buffer, max_y + buffer)

        ax.set_title(title)

        # Calculate the x-tick positions for each chromosome
        x_tick_positions = [(x_positions[i][0] + x_positions[i][regions_number - 1]) / 2 for i in
                            range(num_chromosomes)]
        ax.set_xticks(x_tick_positions)
        ax.set_xticklabels(x_labels, rotation=rotation)

        # Place the legend to the side and divide into columns based on the number of regions
        if legend_labels is None:
            legend_labels = [f'R{i + 1}' for i in range(regions_number)]

        # Adjust the number of columns for the legend
        num_columns = max(1, len(legend_labels) // (X_SIZE // 2))  # Adjust 5 based on the size of the figure

        # Create a custom legend to the side of the plot
        ax.legend(legend_labels, title="Region color", loc='center left', bbox_to_anchor=(1, 0.5), ncol=num_columns)

        if save:
            Graphs._savefig(title, name)
        plt.show()

    """
    This function graphs the CGR representation in 3D to visualize the density of the boxes.
    The nucleotides are in the same position but rotated 45° clockwise.
    The function receives the X, Y starting and ending numbers, the epsilon value (box size)
    and the count_matrix, which is the matrix where [i,j] value is the amount of points in the box [i,j]
    """

    @staticmethod
    def graph_3d_cgr(count_matrix: list[list[int]], name, title='Multifractal measure representation', x_start=0, x_end=1,
                     y_start=0, y_end=1,
                     epsilon=0.01, save=True):
        # set up the figure and axes
        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111, projection='3d')

        _x = np.arange(x_start, x_end, epsilon)
        _y = np.arange(y_start, y_end, epsilon)
        _xx, _yy = np.meshgrid(_x, _y)
        x, y = _xx.ravel(), _yy.ravel()

        # Flatten the count_matrix into a 1D array
        Z = count_matrix.flatten()
        bottom = np.zeros_like(Z)
        width = epsilon
        depth = epsilon

        ax1.bar3d(x, y, bottom, width, depth, Z, shade=True)
        ax1.set_title(title)
        if save:
            Graphs._savefig(title, name)
        plt.show()

    @staticmethod
    def graph_cgr(x_coords, y_coords, name, grid_size=None, epsilon=None, marker='o', markersize=0.005,
                  title="Chaos Game Representation", x_label="X", y_label="Y",
                  x_start=0, x_end=1, y_start=0, y_end=1, with_grid=True, linecolor='red', linewidth=0.1, save=True):
        # Plot the points
        plt.figure(figsize=(10, 10))
        plt.plot(x_coords, y_coords, marker, markersize=markersize)
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        # Set aspect ratio and remove axis
        plt.gca().set_aspect('equal', adjustable='box')
        # plt.axis('off')

        if with_grid:
            square_width = epsilon
            # Plot grid lines for each square
            for i in range(grid_size + 1):
                plt.axvline(i * square_width, color=linecolor, linewidth=linewidth)
                plt.axhline(i * square_width, color=linecolor, linewidth=linewidth)

        if save:
            Graphs._savefig(title, name)

        plt.show()

    @staticmethod
    def graph_coverage_whole(values, title, name, graph_with='seaborn', save=True):
        # Create a list to store colors for each position in the binary sequence
        colors = ['green' if bit == 1 else 'red' for bit in values]

        if graph_with == 'matplotlib':
            plt.figure(figsize=(10, 1))
            plt.bar(range(len(values)), [1] * len(values), color=colors, edgecolor='none', width=1)  # Set width to 1
            plt.title(title)
            plt.xlabel('Position')
            plt.yticks([])
            if save:
                Graphs._savefig(title, name)
            plt.show()

        if graph_with == 'seaborn':
            plt.figure(figsize=(10, 1))
            sns.barplot(x=list(range(len(values))), y=[1] * len(values), palette=colors, edgecolor=None,
                        width=1)  # Use width instead of bar_width
            plt.title(title)
            plt.xlabel('Position')
            plt.ylabel('')
            # Set the number of x-axis ticks
            plt.xticks(
                range(0, len(values), len(values) // 10))  # Adjust the denominator to control the number of ticks
            if save:
                Graphs._savefig(title, name)
            plt.show()

    @staticmethod
    def graph_coverage_without_scale(values: list[int], title, name, save=True):
        total = np.sum(np.abs(values))
        proportions = np.abs(values) / total

        colors = ['green' if value > 0 else 'red' for value in values]

        plt.figure(figsize=(10, 2))
        plt.bar(range(len(values)), np.abs(values), color=colors, width=1.0)

        plt.title(title)
        plt.xlabel('Index')
        plt.yticks([])  # Hide y-axis ticks
        if save:
            Graphs._savefig(title, name)
        plt.show()

    @staticmethod
    def graph_coverage(values: list[int], sequence_name: str, name, save=True):
        total = np.sum(np.abs(values))
        proportions = np.abs(values) / total

        colors = ['green' if value > 0 else 'red' for value in values]

        # Calculate cumulative sum of absolute values
        cumulative_values = np.cumsum(np.abs(values))
        cumulative_values_proportion = cumulative_values / total
        percentage_values = [f"{round(value * 100, 2)}%" for value in cumulative_values_proportion]

        plt.figure(figsize=(15, 10))
        plt.bar(range(len(values)), np.abs(values), color=colors)

        # Set x-axis labels to represent cumulative sum
        plt.xticks(range(len(values)), percentage_values)

        # Display only specific ticks on the x-axis (e.g., every "x" ticks)
        divisions = 4
        xticks_step = len(values) // divisions if len(values) // divisions > 0 else 1
        plt.xticks(range(0, len(values), xticks_step))

        title = f"Coverage of the sequence {sequence_name}"
        plt.title(title)
        plt.xlabel('Representative percentage of the sequence')
        plt.ylabel('Absolute Value')
        if save:
            Graphs._savefig(title, f"{name}/coverage")
        plt.show()

        Graphs.graph_stacked_coverage(values, sequence_name, name, True)

    @staticmethod
    def graph_stacked_coverage(values: list[int], sequence_name: str, name, save=True):
        fig, ax = plt.subplots(figsize=(5, 500))

        # Separate positive and negative values for stacking
        positive_values = [max(0, value) for value in values]
        negative_values = [min(0, value) for value in values]

        # Initialize the starting point of the bars
        current_position = 0

        # Create intercalary sections of green and red stacks
        for positive, negative in zip(positive_values, negative_values):
            ax.bar(0, positive, color='green', bottom=current_position)
            current_position += positive
            ax.bar(0, abs(negative), color='red', bottom=current_position)
            current_position += abs(negative)

        # Set labels and title
        title = f"Stacked coverage of the sequence {sequence_name}"
        ax.set_ylabel(title)
        ax.set_title('Order (bottom-up)')

        if save:
            Graphs._savefig(title, f"{name}/stacked_coverage")
        plt.show()


    @staticmethod
    def graph_vertical_coverage(values: list[int], sequence_name: str, name, save=True):
        total = np.sum(np.abs(values))
        proportions = np.abs(values) / total

        colors = ['green' if value > 0 else 'red' for value in values]

        plt.figure(figsize=(10, 15))
        # plt.barh(range(len(values)), 1, color=colors, height=1, edgecolor='none', align='edge', left=0, linewidth=0, alpha=0.5)
        plt.barh(range(len(values)), proportions, color=colors, height=0.5, edgecolor='none', align='edge', left=0,
                 linewidth=0)

        plt.title('Proportional Color Coding')
        plt.xlabel('Absolute Value')
        plt.ylabel('Index')
        title = f"Coverage of the sequence {sequence_name}"
        plt.title(title)
        if save:
            Graphs._savefig(title, name)
        plt.show()

    @staticmethod
    def graph_linear_fit(fq_values: list[dict], epsilons: list[float], sequence_name: str, name, save=True):
        for index, fq_value in enumerate(fq_values):
            if fq_value['q'] % 2 == 0:
                plt.plot(np.log(epsilons), fq_value['fq'], label=f"q = {fq_value['q']}")
                plt.text(np.log(epsilons[index % len(epsilons)]), fq_value['fq'][index % len(epsilons)],
                         f"q = {fq_value['q']}", va='center', ha='left', fontsize=5)
                plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=2)

        title = f"fq vs ln(ε) for {sequence_name} - {name}"
        plt.title(title)
        if save:
            Graphs._savefig(title, name)
        plt.show()


    @staticmethod
    def graph_linear_regression(fq_values: list[dict], epsilons: list[float], sequence_name: str, name, save=True):
        linear_coefficients = np.polyfit(np.log(epsilons), fq_values['fq'], 1)
        # Create scatter plot
        plt.scatter(np.log(epsilons), fq_values['fq'], label='Data')

        # Generate points for the regression line
        x_values = np.linspace(min(np.log(epsilons)), max(np.log(epsilons)), 100)
        y_values = np.polyval(linear_coefficients, x_values)

        # Plot the regression line
        plt.plot(x_values, y_values, color='red', label='Linear Regression')

        # Add labels and legend
        plt.xlabel('Log(Epsilons)')
        plt.ylabel('fq')
        plt.legend()
        title = f"linear regression for {sequence_name} - {name}"
        if save:
            Graphs._savefig(title, name)

        # Show plot
        plt.show()


    ## Repeats graphs
    # merged with partitions
    @staticmethod
    def _create_partitions(df, size: int, amount_partitions: int, start_col_name: str, length_col_name: str):
        partition_size = size // amount_partitions
        repeat_lengths = np.zeros(amount_partitions)

        for i, row in df.iterrows():
            index = row[start_col_name] // partition_size
            if index >= amount_partitions: #the residue is added to the last partition
                index -= 1
            repeat_lengths[index] += row[length_col_name]
        return repeat_lengths

    @staticmethod
    def graph_distribution_of_repeats_merged(df: pd.DataFrame, size: int, partitions: int = 300,
                                             filter_string: str = None,
                                             filter_column: str = None, regions: int = 3,
                                             plot_type: str = "line", save: bool=True, name: str=None,
                                             filename: str=None):
        if filter_string:
            df = Graphs._filter(df, filter_string, filter_column)

        plt.figure(figsize=(40, 6))

        repeat_lengths = Graphs._create_partitions(df, size, partitions, "query_begin", "repeat_length")

        # Convert repeat_lengths to a numpy array for plotting
        repeat_lengths = np.array(repeat_lengths)

        if plot_type == "line":
            # Plotting as lines
            plt.plot(repeat_lengths, color='gray')
        elif plot_type == "bar":
            # Plotting as bars
            for i, repeat_sum in enumerate(repeat_lengths):
                plt.bar(i, repeat_sum, color='gray')

        title = f"Distribution of Repeats across Sequence {filename} - {name}"
        plt.ylabel("Length of Repeat (bp)")
        plt.xlabel("Repeat")
        plt.title(title)
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        # Add vertical dotted lines at 1/3 and 2/3 of X axis
        x_ticks = np.arange(0, len(repeat_lengths), 1)
        for x in range(1, regions):
            plt.axvline(x=x * max(x_ticks) / regions, color='r', linestyle='--', linewidth=2)

        plt.tight_layout()
        if save:
            Graphs._savefig(title, f"{name}/repeats/RM/distribution_merged")
        plt.show()

    @staticmethod
    def _filter(data, filter_string=None, filter_column="class_family"):
        return data[~data[filter_column].str.contains(filter_string)]

    @staticmethod
    def graph_frequency_of_repeats_grouped(data, col=None, filtering=False, filter_string=None,
                                           filter_column=None, n_max=10, save: bool=True, name: str=None,
                                           filename: str=None):
        grouped_data_sorted = Graphs._group_columns(data, col, "repeat_length", filtering, filter_string, filter_column)
        grouped_data_sorted = grouped_data_sorted.head(n_max)

        title = f"Frequency of Repeats {col} Across Sequence {filename} by {col} - {name}"
        # Plot the distribution of repeats for each class/family
        plt.figure(figsize=(10, 6))
        plt.bar(grouped_data_sorted[col], grouped_data_sorted["repeat_length"])
        plt.xlabel(f"{col}")
        plt.ylabel("Total Repeat Length")
        plt.title(title)
        plt.xticks(rotation=90, fontsize=10)  # Use calculated fontsize
        plt.tight_layout()
        if save:
            Graphs._savefig(title, f"{name}/repeats/RM/frequency_grouped")
        plt.show()

    @staticmethod
    def _group_columns(data, col, sum_col, filtering=False, filter_string=None, filter_column=None):
        if filtering:
            # Filter out rows where the column contains the specified string
            data = Graphs._filter(data, filter_string, filter_column)
        # Group repeats by column and calculate their total length
        grouped_data = data.groupby(col)[sum_col].sum().reset_index()
        grouped_data_sorted = grouped_data.sort_values(by=sum_col, ascending=False)
        return grouped_data_sorted

    @staticmethod
    def graph_distribution_of_repeats(df, col, legend=True, plot_type="line", limit=20,
                                      regions=3, save=True, name=None, filename=None):
        grouped_data_sorted = Graphs._group_columns(df, col, "repeat_length")
        if df.empty or grouped_data_sorted.empty:
            logger.warning("Empty DataFrame provided for Distribution of repeats graph; skipping plot.")
            return

        # Extract unique class/family values
        unique_class_family = grouped_data_sorted.head(limit)[col].tolist()

        # Generate a list of distinct colors dynamically
        num_colors = len(unique_class_family)
        palette = sns.color_palette("hls", num_colors)

        # Dictionary to store class/family and corresponding color
        color_dict = dict(zip(unique_class_family, palette))

        plt.figure(figsize=(30, 6))
        max_value = [0, ""]  # Variable to store the maximum value graphed

        plotted = False
        if plot_type == "line":
            for label in unique_class_family:
                repeat_lengths = np.zeros(len(df))
                for i, row in df.iterrows():
                    if row[col] == label:
                        repeat_lengths[i] = row["repeat_length"]
                        max_value = [max(max_value[0], row["repeat_length"]),
                                     row['repeat'] + " - " + row['class_family']]  # Update the maximum value
                        plotted = True
                plt.plot(repeat_lengths, color=color_dict.get(label), label=label)

        elif plot_type == "bar":
            for i, row in df.iterrows():
                if row[col] in unique_class_family:
                    label = row[col]
                    plt.bar(i, row["repeat_length"], color=color_dict.get(label))
                    max_value = [max(max_value[0], row["repeat_length"]),
                                 row['repeat'] + " - " + row['class_family']]  # Update the maximum value
                    plotted = True

        title = f"Distribution of Repeats Across Sequence {filename} by {col} - {name}"
        plt.ylabel("Length of Repeat (bp)")
        plt.xlabel("Repeat")
        plt.title(title)
        logger.info(f"Distribution of repeats - max value {max_value}")
        if max_value[0] > 0:
            plt.ylim(0, max_value[0] * 1.1)
        else:
            logger.warning("No valid data to plot; skipping ylim setting.")

        x_ticks = np.arange(0, len(df), 1)
        for x in range(1, regions):
            plt.axvline(x=x * max(x_ticks) / regions, color='r', linestyle='--', linewidth=2)

        # Create legend using unique class/family names
        if legend and plotted:
            plt.legend(ncol=2, loc='upper left', bbox_to_anchor=(1, 1))
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        if save:
            Graphs._savefig(title, f"{name}/repeats/RM/distribution_separately")
        plt.show()

    @staticmethod
    def graph_distribution_of_repeats_subplots(df, col="class_family", legend=True, limit=20,
                                               regions=3, shared_y_axis=False, save=True, name=None,
                                               filename=None):
        grouped_data_sorted = Graphs._group_columns(df, col, "repeat_length")
        # Extract unique class/family values
        unique_class_family = grouped_data_sorted.head(limit)[col].tolist()

        # Generate a list of distinct colors dynamically
        num_colors = len(unique_class_family)
        palette = sns.color_palette("hls", num_colors)

        # Dictionary to store class/family and corresponding color
        color_dict = dict(zip(unique_class_family, palette))

        num_subplots = len(unique_class_family)
        fig, axs = plt.subplots(num_subplots, 1, figsize=(30, 3 * num_subplots))

        max_y = 0  # Variable to store the maximum y-axis value among all subplots

        for i, label in enumerate(unique_class_family):
            ax = axs[i]
            repeat_lengths = np.zeros(len(df))
            for i, row in df.iterrows():
                if row[col] == label:
                    repeat_lengths[i] = row["repeat_length"]
            ax.plot(repeat_lengths, color=color_dict.get(label), label=label)
            ax.set_ylabel("Length of Repeat (bp)")
            ax.set_xlabel("Repeat")
            ax.set_title(f"Distribution of {label} repeats across Sequence {filename} - {name}")
            ax.grid(axis='y', linestyle='--', alpha=0.7)
            ax.set_ylim(0, np.max(repeat_lengths) * 1.1)  # Set the y-axis limit to 110% of the maximum value graphed

            x_ticks = np.arange(0, len(df), 1)
            for x in range(1, regions):
                ax.axvline(x=x * max(x_ticks) / regions, color='r', linestyle='--', linewidth=2)

            if legend:
                ax.legend()

            if shared_y_axis:
                max_y = max(max_y, np.max(repeat_lengths))  # Update the maximum y-axis value

        if shared_y_axis:
            # Set the same y-axis limits for all subplots
            for ax in axs:
                ax.set_ylim(0, max_y)

        plt.tight_layout()
        title = f"Distribution of repeats across Sequence {filename} by {col}"
        if save:
            Graphs._savefig(title, f"{name}/repeats/RM/distribution_subplots")
        plt.show()


    @staticmethod
    def graph_recursive_repeats_largest_values_from_database(df, col=None, n_max=None, save=True, name=None,
                                                filename=None):
        df = df.sort_values(by='largest_value', ascending=False)
        if n_max is not None:
            df = df.head(n_max)

        plt.figure(figsize=(10, 6))
        plt.bar(df[col], df['largest_value'], color='blue')

        title = f"Recursively found repeats - largest values in Sequence {filename} - {name}"
        plt.xlabel('Sequence')
        plt.ylabel('Largest Value')
        plt.title(title)
        plt.xticks(rotation=90)

        plt.tight_layout()
        if save:
            Graphs._savefig(title, f"{name}/repeats/recursive/largest_values")
        plt.show()

    @staticmethod
    def graph_grouped_by_recursive_repeat_length(df, col=None, save=True, name=None, filename=None):
        grouped = df.groupby('repeat_length')

        plt.figure(figsize=(14, 8))

        for repeat_length, group in grouped:
            group = group.sort_values(by='largest_value', ascending=False)
            plt.bar(group[col], group['largest_value'], label=f'Repeat Length: {repeat_length}')

        title = f"Largest Values Grouped by Repeat Length of Sequence {filename} - {name}"
        plt.xlabel('Sequence')
        plt.ylabel('Largest Value')
        plt.title(title)
        plt.xticks(rotation=90)
        plt.legend(title='Repeat Length')

        plt.tight_layout()
        if save:
            Graphs._savefig(title, f"{name}/repeats/recursive/by_repeat_length")
        plt.show()

    @staticmethod
    def graph_individual_plots_by_recursive_repeat_length(df, col=None, save=True, name=None,
                                                          filename=None):
        grouped = df.groupby('repeat_length')

        for repeat_length, group in grouped:
            group_sorted = group.sort_values(by='largest_value', ascending=False)

            plt.figure(figsize=(10, 6))
            plt.bar(group_sorted[col], group_sorted['largest_value'], color='blue')

            title = f"Largest Values for Repeat Length {repeat_length} of Sequence {filename} - {name}"
            plt.xlabel('Sequence')
            plt.ylabel('Largest Value')
            plt.title(title)
            plt.xticks(rotation=90)

            plt.tight_layout()

            if save:
                route = f"{name}/repeats/recursive/repeat_length_{repeat_length}"
                Graphs._savefig(title, route)
            plt.show()

    @staticmethod
    def graph_distribution_of_genes_merged(df, name: str, size: int, partitions: int,regions: int, plot_type: str,
                                           chromosome_name: str, save: bool):
        df = df[df['feature'] == 'gene']
        df = df.reset_index(drop=True)
        plt.figure(figsize=(40, 6))
        lengths = Graphs._create_partitions(df, size, partitions, start_col_name="start_position",
                                            length_col_name="length")

        # Convert lengths to a numpy array for plotting
        lengths = np.array(lengths)

        if plot_type == "line":
            plt.plot(lengths, color='gray')
        elif plot_type == "bar":
            for i, repeat_sum in enumerate(lengths):
                plt.bar(i, repeat_sum, color='gray')

        title = f"Distribution of Genes Across Sequence {chromosome_name} - {name}"
        plt.ylabel("Length of Gene (bp)")
        plt.xlabel("Genes")
        plt.title(title)
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        # Add vertical dotted lines at 1/3 and 2/3 of X axis
        x_ticks = np.arange(0, len(lengths), 1)
        for x in range(1, regions):
            plt.axvline(x=x * max(x_ticks) / regions, color='r', linestyle='--', linewidth=2)

        plt.tight_layout()
        if save:
            route = f"{name}/genes/gtf_merged"
            Graphs._savefig(title, route)
        plt.show()


    @staticmethod
    def graph_distribution_of_genes(df, name: str, legend: bool, plot_type: str, limit: int, regions: int,
                                    chromosome_name: str, save: bool):

        grouped_data_sorted = Graphs._group_columns(df, "feature", "length")

        unique_type = grouped_data_sorted.head(limit)["feature"].tolist()

        num_colors = len(unique_type)
        palette = sns.color_palette("hls", num_colors)

        color_dict = dict(zip(unique_type, palette))

        plt.figure(figsize=(30, 6))
        max_value = [0, ""]  # Variable to store the maximum value graphed

        if plot_type == "line":
            for label in unique_type:
                repeat_lengths = np.zeros(len(df))
                for i, row in df.iterrows():
                    if row["feature"] == label:
                        repeat_lengths[i] = row["length"]
                        max_value = [max(max_value[0], row["length"]), row['feature']]  # Update the maximum value
                plt.plot(repeat_lengths, color=color_dict.get(label), label=label)

        elif plot_type == "bar":
            for i, row in df.iterrows():
                if row["feature"] in unique_type:
                    label = row["feature"]
                    plt.bar(i, row["length"], color=color_dict.get(label))
                    max_value = [max(max_value[0], row["length"]), row['feature']]  # Update the maximum value

        title = f"Distribution of Genes Across Sequence {chromosome_name} - {name}"
        plt.ylabel("Length of Gene (bp)")
        plt.xlabel("Genes")
        plt.title(title)
        plt.ylim(0, max_value[0] * 1.1)  # Set the y-axis limit to 110% of the maximum value graphed

        x_ticks = np.arange(0, len(df), 1)
        for x in range(1, regions):
            plt.axvline(x=x * max(x_ticks) / regions, color='r', linestyle='--', linewidth=2)

        # Create legend using unique class/family names
        if legend:
            plt.legend(ncol=2, loc='upper left', bbox_to_anchor=(1, 1))
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        if save:
            route = f"{name}/genes/gtf"
            Graphs._savefig(title, route)
        plt.show()
