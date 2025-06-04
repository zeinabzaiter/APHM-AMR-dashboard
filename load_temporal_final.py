import os
import pandas as pd
import sqlalchemy

def load_temporal_final(sqlite_path="db/amr_data.sqlite", data_folder="data"):
    """
    Parcourt les fichiers sources (Staph aureus et Enterococcus faecium),
    calcule le %R hebdomadaire, la moyenne mobile 8 semaines, l'écart-type,
    les intervalles de confiance (±1 écart-type) et détecte les outliers.
    Insère le résultat dans la table 'resistance_temps' de la base SQLite.
    """
    engine = sqlalchemy.create_engine(f"sqlite:///{sqlite_path}")
    all_rows = []

    # --- Chargement Staphylococcus aureus ---
    staph_path = os.path.join(data_folder, "Export_StaphAureus_COMPLET.csv")
    if os.path.exists(staph_path):
        df_staph = pd.read_csv(staph_path, sep=",")
        # Renommer colonnes
        rename_map = {
            'Week': 'week_number',
            'Vancomycine': 'Vancomycin',
            'Teicoplanine': 'Teicoplanin',
            'Gentamycine': 'Gentamicin',
            'Oxacilline': 'Oxacillin',
            'Daptomycine': 'Daptomycin',
            'Dalbavancine': 'Dalbavancin',
            'Clindamycine': 'Clindamycin',
            'Cotrimoxazole': 'Cotrimoxazole',
            'Linezolide': 'Linezolid'
        }
        df_staph = df_staph.rename(columns=rename_map)
        # Date ISO pour 2024
        df_staph['week'] = pd.to_datetime('2024-' + df_staph['week_number'].astype(str) + '-1', format='%G-%V-%u')
        # Format long
        antibiotiques = ['Vancomycin', 'Teicoplanin', 'Gentamicin', 'Oxacillin',
                         'Daptomycin', 'Dalbavancin', 'Clindamycin', 'Cotrimoxazole', 'Linezolid']
        staph_long = df_staph.melt(
            id_vars=['week', 'lib_germe'],
            value_vars=antibiotiques,
            var_name='antibiotic',
            value_name='susceptibility'
        )
        staph_long = staph_long[staph_long['lib_germe'] == 'Staphylococcus aureus']
        # Calcul % résistants
        percent_staph = (
            staph_long
            .dropna(subset=['susceptibility'])
            .assign(resistant=lambda df: df['susceptibility'].str.upper() == 'R')
            .groupby(['week', 'antibiotic'], as_index=False)
            .agg(total_isolates=('resistant', 'size'), total_resistants=('resistant', 'sum'))
        )
        percent_staph['percent_resistant'] = 100 * percent_staph['total_resistants'] / percent_staph['total_isolates']
        percent_staph['bacterie'] = 'Staphylococcus aureus'
        all_rows.append(percent_staph)

    # --- Chargement Enterococcus faecium ---
    entero_path = os.path.join(data_folder, "Enterococcus_faecium_groupes_antibiotiques.xlsx")
    if os.path.exists(entero_path):
        df_entero = pd.read_excel(entero_path, sheet_name=0)
        # Renommage colonnes
        rename_map2 = {
            'Numéro semaine': 'week_number',
            'Ampicilline': 'Ampicillin',
            'Vancomycine': 'Vancomycin',
            'Teicoplanine': 'Teicoplanin',
            'Gentamicine': 'Gentamicin',
            'Linezolide': 'Linezolid',
            'Daptomycine': 'Daptomycin',
            'Tigecycline': 'Tigecycline'
        }
        df_entero = df_entero.rename(columns=rename_map2)
        df_entero['week'] = pd.to_datetime('2024-' + df_entero['week_number'].astype(str) + '-1', format='%G-%V-%u')
        antibiotiques_entero = ['Ampicillin', 'Vancomycin', 'Teicoplanin',
                                'Gentamicin', 'Linezolid', 'Daptomycin', 'Tigecycline']
        entero_long = df_entero.melt(
            id_vars=['week', 'LIB_GERME'],
            value_vars=antibiotiques_entero,
            var_name='antibiotic',
            value_name='susceptibility'
        )
        entero_long = entero_long[entero_long['LIB_GERME'] == 'Enterococcus faecium']
        percent_entero = (
            entero_long
            .dropna(subset=['susceptibility'])
            .assign(resistant=lambda df: df['susceptibility'].str.upper() == 'R')
            .groupby(['week', 'antibiotic'], as_index=False)
            .agg(total_isolates=('resistant', 'size'), total_resistants=('resistant', 'sum'))
        )
        percent_entero['percent_resistant'] = 100 * percent_entero['total_resistants'] / percent_entero['total_isolates']
        percent_entero['bacterie'] = 'Enterococcus faecium'
        all_rows.append(percent_entero)

    # Concaténer
    if all_rows:
        df_all = pd.concat(all_rows, ignore_index=True)
        df_all = df_all.sort_values(['bacterie', 'antibiotic', 'week']).reset_index(drop=True)

        # Calcul moyenne mobile 8 semaines, écart-type et IC
        df_all['moyenne_mobile'] = df_all.groupby(['bacterie', 'antibiotic'])['percent_resistant']             .transform(lambda x: x.rolling(window=8, min_periods=1).mean())
        df_all['ecart_type'] = df_all.groupby(['bacterie', 'antibiotic'])['percent_resistant']             .transform(lambda x: x.rolling(window=8, min_periods=1).std().fillna(0))

        # IC supérieur / inférieur = moyenne ± 1 écart-type
        df_all['ic_sup'] = df_all['moyenne_mobile'] + df_all['ecart_type']
        df_all['ic_inf'] = df_all['moyenne_mobile'] - df_all['ecart_type']

        # Détecter outliers : %R en dehors de l'IC
        df_all['outlier'] = (df_all['percent_resistant'] > df_all['ic_sup']) |                             (df_all['percent_resistant'] < df_all['ic_inf'])

        # Choisir colonnes à insérer
        df_final = df_all[[
            'bacterie', 'antibiotic', 'week', 'percent_resistant',
            'moyenne_mobile', 'ecart_type', 'ic_sup', 'ic_inf', 'outlier'
        ]]

        # Charger dans SQLite
        with engine.connect() as conn:
            df_final.to_sql("resistance_temps", conn, if_exists="replace", index=False)
        print("Table 'resistance_temps' mise à jour dans la base SQLite (2024).")
    else:
        print("Aucun fichier source valide trouvé dans 'data/'.")

if __name__ == "__main__":
    load_temporal_final()
