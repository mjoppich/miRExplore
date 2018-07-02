from collections import defaultdict

import os


class feedbackDB:


    def __init__(self, filepath):

        self.all_feedback = []
        self.filepath = filepath

        self.feedbackByID = defaultdict(list)

        if self.filepath != None:

            if os.path.isfile(filepath):
                with open(filepath, 'r') as fin:

                    for line in fin:
                        line = line.strip()
                        aline = eval(line)
                        self.all_feedback.append(aline)

                        self.feedbackByID[aline[1]].append(aline)


    def get_feedback(self, dataid):
        return self.feedbackByID.get(dataid, None)

    def add_feedback(self, feedback):
        self.all_feedback.append(feedback)
        self.feedbackByID[feedback[1]].append(feedback)

        self.save_to_file()

    def save_to_file(self):

        if self.filepath != None:
            with open(self.filepath, 'w') as fout:
                for elem in self.all_feedback:
                    fout.write(str(elem) + "\n")

