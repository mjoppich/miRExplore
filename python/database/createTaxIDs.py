from database.Neo4JInterface import neo4jInterface

db = neo4jInterface(simulate=False)
db.createUniquenessConstraint('TAX', 'id')


db.createNodeIfNotExists(["TAX"], {'id': 9606, 'name':'Homo sapiens', 'short':'hsa'})
db.createNodeIfNotExists(["TAX"], {'id': 10116, 'name':'Mus musculus', 'short':'mmu'})

db.close()