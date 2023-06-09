import mysql.connector
import configparser

## Abandoned project ignore, but working!
# Mysql connection use case, no specific modules


class DB:
    def __init__(self, profile: str, cfg: str) -> None:
        config = configparser.ConfigParser()
        config.read(cfg)
        self.database = config[profile]['database']
        self.user = config[profile]['user']
        self.password = config[profile]['password']
        self.host = config[profile]['host']
        self.cnx = None
        self.cursor = None
        self.current_table = None

    def connect(self) -> None:
        self.cnx = mysql.connector.connect(user=self.user, password=self.password, host=self.host, database=self.database)
        self.cursor = self.cnx.cursor()

    def close(self) -> None:
        self.cnx.close()

    def cursor(self):
        return self.cursor

    def set_current_table(self, table):
        self.current_table = table

    def add_new_columns(self, new_columns) -> None:
        result = f'ALTER TABLE `{self.current_table}` ADD '
        for column in new_columns:
            result += f'`{column}` VARCHAR(254)'
        print(f''' ALTER TABLE `Feature` ADD `{new_columns[0]}` VARCHAR(254) NOT NULL AFTER `end`, ADD `{new_columns[1]}` VARCHAR(254) NOT NULL AFTER `{new_columns[0]}`;''')



db = DB('sRNA', "params.ini")
db.connect()
newColumns = ["sequence_assembly_key", "sequence_assembly_value"]
db.add_new_columns(newColumns)



target = "Genome"
origin = "Feature"
database = "sRNAPlantPortal"

targetPK = ""



query = f'''SHOW INDEX FROM {database}.{target} WHERE Key_name = 'PRIMARY';'''

cursor = db.cursor()
db.cursor.execute(query)
for line in cursor:
    targetPK = line[4]

print(query)

# target extract pk + new pKs

# origin add new Fks

# match

# insert data in columns

db.close()



