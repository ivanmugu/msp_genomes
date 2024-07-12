import os

class Tidy:
    def __init__(self, fasta_file, force_overwrite = False):
        self.fasta_file = fasta_file
        self.force_overwrite = force_overwrite
                
    @staticmethod
    def get_lines(fasta_file:str) -> list[str]:
        """Extracts the lines in the fasta file

        Args:
            fasta_file (str): Path to fasta file

        Returns:
            list[str]: List with lines
        """
        with open(fasta_file, 'r') as f:
            lines = f.readlines()
        return lines
    
    def overwrite(self, path:str ) -> str:
        """If force_overwrite is True, the path is returned as is. Otherwise, the path is modified to avoid overwriting the original file.

        Args:
            path (str): _description_

        Returns:
            str: Path to save.
        """
        if self.force_overwrite:
            modified_path = path
        else:
            splits_base = os.path.basename(path).split(".")
            modified_base = splits_base[0] + "_copy"
            modified_path = os.path.join(os.path.dirname(path), modified_base + "." + ".".join(splits_base[1:]))
        return modified_path

    def join_header_by_underscore(self) -> None:
        """Refers to Issue 9 on virulencefinder_db repo (https://bitbucket.org/genomicepidemiology/virulencefinder_db/issues/9/keyerror-when-running-virulencefinder)\n
            Replace spaces in headers with underscores in fasta file"
        """
        lines = self.get_lines(self.fasta_file)
        write_path = self.overwrite(self.fasta_file)
        with open(write_path, 'w') as f:
            for line in lines:
                if line.startswith('>'):
                    f.write(line.replace(' ', '_'))
                else:
                    f.write(line)
                    
#Tidy("/home/people/s220672/databases/virulencefinder_db/virulence_ecoli.fsa", force_overwrite = True).join_header_by_underscore()
