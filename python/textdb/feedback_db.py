from collections import defaultdict

import os


class feedbackDB:


    def __init__(self, filepath):

        self.all_feedback = []
        self.filepath = filepath

        if os.path.isfile(filepath):
            with open(filepath, 'r') as fin:

                for line in fin:
                    line = line.strip()
                    aline = eval(line)
                    self.all_feedback.append(aline)

    def add_feedback(self, feedback):
        self.all_feedback.append(feedback)
        self.save_to_file()

    def save_to_file(self):

        with open(self.filepath, 'w') as fout:
            for elem in self.all_feedback:
                fout.write(str(elem) + "\n")

