import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle

import os


class Graphs:

    @staticmethod
    def _savefig(title, name):
        directory = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), "out/graphs")
        os.makedirs(directory, exist_ok=True)

        actual_path = f"{directory}/{name}"
        plt.tight_layout()
        if not os.path.exists(actual_path):
            os.makedirs(actual_path)
        plt.savefig(f'{actual_path}/{title}.png')

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
    def graph_bars_grouped(x_array, y_array, title, name, regions_number=3, y_label=None, x_labels=None, legend_labels=None,
                           regions_colors=None, rotation=45, y_range: list[int] = None, top_labels=False, save=True):
        fig, ax = plt.subplots()
        bar_width = 0.1  # Width of each bar
        num_chromosomes = len(x_array) // regions_number

        if x_labels is None:
            x_labels = [f'Chromosome {i + 1}' for i in range(num_chromosomes)]

        if regions_colors is None:
            colors = ['b', 'r', 'g', 'c', 'y', 'm', 'k', 'w']
            colors = colors[:regions_number]
            colors_cycle = cycle(colors)
        else:
            colors_cycle = cycle(regions_colors)

        # Calculate the x-positions for the bars
        x_positions = [i * (regions_number + 2) + np.arange(regions_number) for i in range(num_chromosomes)]

        for i in range(num_chromosomes):
            for j in range(regions_number):
                color = next(colors_cycle)
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

        # Create a custom legend to clarify which region corresponds to which color
        ax.legend(legend_labels, title="Region color")

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
                Graphs._savefig(title, )
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
        title = f"Coverage of the sequence {sequence_name}"
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

        title = f"fq vs ln(ε) for {sequence_name}"
        plt.title(title)
        if save:
            Graphs._savefig(title, name)
        plt.show()
