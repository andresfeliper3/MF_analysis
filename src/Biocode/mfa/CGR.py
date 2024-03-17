import numpy as np
from src.Biocode.sequences.Sequence import Sequence

from src.Biocode.graphs.Graphs import Graphs

from utils.logger import logger

class CGR:
    def __init__(self, sequence: Sequence, corners=None):
        self.mi_grid = None
        self.sum_nucleotides = None
        self.grid_size = None
        if corners is None:
            corners = {'A': (0, 0), 'C': (0, 1), 'G': (1, 1), 'T': (1, 0)}
        self.sequence = sequence
        self.seq = sequence.get_sequence()
        self.len_seq = len(self.seq)
        self.corners = corners
        self.cover = [0]
        self.cover_percentage = 0.0

    def generate_only_cgr(self, graph=True, title="Chaos Game Representation") -> list[list[float]]:
        # Initialize the current position
        current_position = np.array([0.5, 0.5])

        # Lists to store the x and y coordinates for plotting
        self.x_coords = [current_position[0]]
        self.y_coords = [current_position[1]]

        # Iterate through the sequence and plot points
        last_exists = True if self.seq[0] in self.corners else False
        for nucleotide in self.seq:
            if nucleotide in self.corners:
                corner = self.corners[nucleotide]
                if last_exists:
                    self.cover[-1] += 1
                else:
                    self.cover.append(1)
                last_exists = True
            else:
                if last_exists:
                    self.cover.append(-1)
                else:
                    self.cover[-1] -= 1
                last_exists = False
                continue
            # Calculate the midpoint between the current position and the selected corner
            current_position = (current_position + np.array(corner)) / 2
            self.x_coords.append(current_position[0])
            self.y_coords.append(current_position[1])

        self.sum_nucleotides = sum(x for x in self.cover if x > 0)
        self.cover_percentage = self.sum_nucleotides / self.len_seq
        self.sequence.set_cover(self.cover)
        self.sequence.set_cover_percentage(self.cover_percentage)

        if graph:
            Graphs.graph_cgr(self.x_coords, self.y_coords, title=title, name=self.sequence.get_name(),with_grid=False)

        self.sequence.set_cover(self.cover)
        return [self.x_coords, self.y_coords]

    def generate_cgr_counting_grid_cells(self, graph=True, epsilon=0.01, title="Chaos Game Representation") -> list[
        list[float]]:

        # Initialize the current position
        current_position = np.array([0.5, 0.5])

        # Lists to store the x and y coordinates for plotting
        self.x_coords = [current_position[0]]
        self.y_coords = [current_position[1]]

        # Grid to store Mi values
        self.grid_size = int(1 / epsilon)
        self.mi_grid = np.zeros((self.grid_size, self.grid_size), dtype=int)

        # Iterate through the sequence and plot points
        last_exists = True if self.seq[0] in self.corners else False
        for seq_index, nucleotide in enumerate(self.seq):
            if nucleotide in self.corners:
                corner = self.corners[nucleotide]
                if last_exists:
                    self.cover[-1] += 1
                else:
                    self.cover.append(1)
                last_exists = True
            else:
                if last_exists:
                    self.cover.append(-1)
                else:
                    self.cover[-1] -= 1
                last_exists = False
                continue
            # Calculate the midpoint between the current position and the selected corner
            current_position = (current_position + np.array(corner)) / 2
            self.x_coords.append(current_position[0])
            self.y_coords.append(current_position[1])

            # Calculate the row and column index of the current point
            row = int(current_position[1] // epsilon) if int(current_position[1] // epsilon) < self.grid_size else int(
                current_position[1] // epsilon) - 1
            col = int(current_position[0] // epsilon) if int(current_position[0] // epsilon) < self.grid_size else int(
                current_position[0] // epsilon) - 1

            if seq_index % 1000000 == 0:
                logger.info(f"In nucleotide: {seq_index}")

            # Increment the count for the current grid cell
            try:
                self.mi_grid[row, col] += 1
            except:
                print("Error here")
                print("position:", current_position[0], ";", current_position[1], "epsilon:", epsilon, "row:", row,
                      "col:", col)

        self.sum_nucleotides = sum(x for x in self.cover if x > 0)
        self.cover_percentage = self.sum_nucleotides / self.len_seq
        self.sequence.set_cover(self.cover)
        self.sequence.set_cover_percentage(self.cover_percentage)

        if graph:
            Graphs.graph_cgr(self.x_coords, self.y_coords, self.grid_size, title=title, epsilon=epsilon)

        self.sequence.set_cover(self.cover)
        return self.mi_grid

    def graph_cgr(self, title, name):
        Graphs.graph_cgr(self.x_coords, self.y_coords, name=name, grid_size=self.grid_size, title=title, with_grid=False)

    def get_corners(self):
        return self.corners

    def get_grid_size(self) -> int:
        return self.grid_size

    def get_mi_grid(self) -> list[list[float]]:
        return self.mi_grid

    def get_cover(self) -> list[int]:
        return self.cover

    def get_cover_percentage(self) -> float:
        return self.cover_percentage
