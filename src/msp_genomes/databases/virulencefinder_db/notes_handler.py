import os

dir = os.path.dirname(__file__)
with open(os.path.join(dir, "notes.txt")) as f:
    contents = f.read()
    
contents = contents.replace("\t", " ")

with open('notes.txt', 'w') as file:
    file.write(contents)
    

