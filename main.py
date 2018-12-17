import EasySectionAnalysis as section

def shit():
    print("Shit")

def main():
    shit()
    wingbox = section.WinboxSection()
    wingbox.createGeometry()
    wingbox.plotGeometry()
    
if __name__ == "__main__":
    main()