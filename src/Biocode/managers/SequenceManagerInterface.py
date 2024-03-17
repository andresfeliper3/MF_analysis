class SequenceManagerInterface:
    def generate_mfa_generator(self):
        """Generate the MFA class"""
        pass

    def generate_mfa(self):
        """Generate mfa_results with the MFA values"""
        pass

    def generate_degree_of_multifractality(self):
        """Generate degrees of multifractality"""
        pass

    def graph_cgr(self):
        "Generate CGR without the grid"
        pass

    def graph_3d_cgr(self, grid_size):
        """Graph CGR density in a 3D chart"""
        pass

    def graph_linear_fit(self):
        """Graph linear fit for fq vs ln(epsilon)"""
        pass

    def graph_multifractal_spectrum(self):
        """Graph multifractal spectrum of Dq vs q"""
        pass

    def graph_correlation_exponent(self):
        """Graph t(q) vs q"""
        pass

    def _attach_cover_data(self):
        """Attach to the fields of the class the cover and cover_percentage of the sequence"""
        pass

    def calculate_multifractal_analysis_values(self):
        """Generate mfa generators, generate mfa values, attach the degrees of multifractality, the cover and cover
        percentage"""
        self.generate_mfa()
        self.generate_degree_of_multifractality()
        self._attach_cover_data()

    def graph_multifractal_analysis(self, _3d_cgr=False, linear_fit=True, degrees_of_multifractality=None,
                                    multifractal_spectrum=True, correlation_exponent=True):
        """Graph the multifractal analysis values using the chart of 3D density of points, the linear fit fq vs q,
        the multifractal spectrum Dq vs q, and the correlation exponent t(q) vs q"""
        if _3d_cgr:
            self.graph_3d_cgr()
        if linear_fit:
            self.graph_linear_fit()
        if multifractal_spectrum:
            self.graph_multifractal_spectrum()
        if correlation_exponent:
            self.graph_correlation_exponent()

    def graph_coverage(self):
        """Graph the barplot representing the coverage of the DNA sequence (the nucleotides that have been
        identified)"""
        pass

    def calculate_and_graph(self):
        """Generate MFA values and graph them along with the coverage"""
        self.calculate_multifractal_analysis_values()
        self.graph_multifractal_analysis()
        self.graph_coverage()

    def set_cover(self, cover):
        pass

    def set_cover_percentage(self, cover_percentage):
        pass
