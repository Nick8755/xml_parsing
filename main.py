#create xml file
import xml
import xml.etree.ElementTree as ET

tree = ET.parse('full database.xml')
root = tree.getroot()

print(root.tag)
print(root.attrib)

for child in root:
    print(child.tag, child.attrib)