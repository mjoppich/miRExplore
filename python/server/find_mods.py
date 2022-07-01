from modulefinder import ModuleFinder
import sys

finder = ModuleFinder()
finder.run_script(sys.argv[1])

allmods = set()

for name, mod in finder.modules.items():
        name = name.split(".")[0]
        
        if name[0] == '_':
            continue

        allmods.add(name)

for x in sorted(allmods):
    print(x)
