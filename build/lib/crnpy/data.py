# crnpy/data.py
"""## crnpy.data
Data module for crnpy.

This module contains data for the crnpy package.

Attributes:
    cutoff_rigidity (list): Cutoff rigidity values for the whole world. See [crnpy.crnpy.cutoff_rigidity][]
    neutron_detectors (list): Neutron detector locations. See [crnpy.crnpy.find_neutron_monitor][]

"""

cutoff_rigidity = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                   [0.0, 0.0, 0.02, 0.02, 0.02, 0.03, 0.03, 0.04, 0.03, 0.03, 0.04, 0.03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                   [0.04, 0.07, 0.09, 0.08, 0.13, 0.08, 0.16, 0.14, 0.13, 0.15, 0.17, 0.13, 0.08, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02, 0.05, 0.04],
                   [0.18, 0.28, 0.34, 0.32, 0.35, 0.41, 0.38, 0.42, 0.47, 0.47, 0.46, 0.43, 0.33, 0.21, 0.11, 0.04, 0.0, 0.0, 0.0, 0.0, 0.02, 0.08, 0.15, 0.2, 0.18],
                   [0.6, 0.66, 0.74, 0.75, 0.76, 0.79, 0.81, 0.87, 0.91, 0.98, 1.01, 0.96, 0.72, 0.55, 0.3, 0.19, 0.07, 0.04, 0.05, 0.06, 0.15, 0.31, 0.45, 0.58, 0.6],
                   [1.13, 1.3, 1.35, 1.39, 1.39, 1.44, 1.47, 1.56, 1.72, 1.76, 1.86, 1.8, 1.42, 1.08, 0.74, 0.48, 0.32, 0.21, 0.22, 0.31, 0.49, 0.7, 0.94, 1.14, 1.13],
                   [2.1, 2.18, 2.22, 2.3, 2.3, 2.38, 2.45, 2.59, 2.75, 2.96, 3.1, 2.97, 2.31, 1.85, 1.39, 0.93, 0.65, 0.5, 0.5, 0.65, 0.94, 1.38, 1.75, 2.02, 2.1],
                   [3.4, 3.49, 3.49, 3.53, 3.58, 3.69, 3.81, 4.05, 4.3, 4.61, 4.72, 4.54, 3.59, 2.91, 2.27, 1.6, 1.18, 0.93, 0.9, 1.21, 1.71, 2.4, 2.96, 3.24, 3.4],
                   [4.98, 5.09, 5.07, 5.12, 5.25, 5.39, 5.53, 5.77, 6.03, 6.3, 6.35, 5.96, 4.96, 4.31, 3.42, 2.68, 1.96, 1.57, 1.56, 1.94, 2.85, 3.92, 4.6, 4.97, 4.98],
                   [7.33, 7.26, 7.11, 7.08, 7.32, 7.67, 7.99, 8.34, 8.93, 9.3, 9.22, 8.63, 6.72, 5.65, 4.88, 4.06, 2.99, 2.46, 2.43, 3.0, 4.24, 5.47, 6.56, 7.24, 7.33],
                   [9.87, 9.73, 9.83, 10.05, 10.51, 10.87, 11.13, 11.02, 11.29, 11.51, 11.16, 10.3, 9.33, 8.11, 6.43, 5.3, 4.31, 3.54, 3.5, 4.29, 5.87, 8.43, 9.54, 9.86, 9.87],
                   [11.73, 11.77, 11.79, 12.13, 12.61, 13.1, 13.58, 13.84, 13.86, 13.66, 13.19, 12.51, 10.83, 9.78, 8.85, 7.24, 5.65, 4.51, 4.48, 5.76, 8.89, 10.72, 11.4, 11.71, 11.73],
                   [13.45, 13.68, 13.9, 14.17, 14.62, 15.05, 15.32, 15.35, 15.19, 14.81, 14.22, 13.52, 12.32, 11.61, 10.73, 9.53, 7.8, 6.22, 6.14, 8.09, 10.88, 12.19, 12.93, 13.35, 13.45],
                   [14.31, 14.63, 14.9, 15.26, 15.76, 16.23, 16.49, 16.44, 16.15, 15.65, 14.99, 14.29, 13.2, 12.61, 11.9, 10.93, 8.99, 6.99, 7.04, 9.87, 12.07, 13.07, 13.73, 14.18, 14.31],
                   [14.71, 15.11, 15.47, 15.93, 16.49, 17.0, 17.25, 17.15, 16.78, 16.21, 15.54, 14.88, 13.9, 13.38, 12.79, 11.96, 10.64, 9.38, 9.64, 11.59, 12.74, 13.54, 14.1, 14.57, 14.71],
                   [14.7, 15.16, 15.62, 16.18, 16.82, 17.36, 17.62, 17.5, 17.08, 16.5, 15.87, 15.29, 14.43, 13.95, 13.44, 12.83, 11.92, 11.1, 11.33, 12.25, 13.07, 13.66, 14.1, 14.54, 14.7],
                   [14.3, 14.8, 15.36, 16.03, 16.75, 17.33, 17.6, 17.48, 17.07, 16.52, 15.97, 15.49, 14.77, 14.34, 13.89, 13.38, 12.76, 12.16, 12.12, 12.57, 13.13, 13.48, 13.75, 14.14, 14.3],
                   [13.56, 14.08, 14.72, 15.5, 16.29, 16.9, 17.18, 17.1, 16.73, 16.24, 15.81, 15.46, 14.9, 14.53, 14.12, 13.68, 13.17, 12.66, 12.43, 12.65, 12.95, 13.04, 13.11, 13.41, 13.56],
                   [12.55, 13.06, 13.78, 14.65, 15.48, 16.08, 16.38, 16.35, 16.06, 15.66, 15.35, 15.15, 14.78, 14.48, 14.14, 13.75, 13.31, 12.83, 12.52, 12.53, 12.57, 12.4, 12.26, 12.43, 12.55],
                   [11.33, 11.84, 12.6, 13.51, 14.33, 14.89, 15.18, 15.23, 15.02, 14.71, 14.54, 14.5, 14.39, 14.2, 13.95, 13.63, 13.25, 12.79, 12.41, 12.24, 12.02, 11.56, 11.15, 11.23, 11.33],
                   [9.92, 10.36, 11.19, 12.09, 12.8, 13.32, 13.56, 13.67, 13.54, 13.32, 13.29, 13.42, 13.69, 13.66, 13.54, 13.32, 13.0, 12.58, 12.14, 11.79, 11.34, 10.67, 9.87, 9.8, 9.92],
                   [8.19, 8.64, 9.37, 10.07, 10.6, 10.89, 11.19, 11.36, 11.19, 10.53, 10.67, 11.73, 12.6, 12.84, 12.91, 12.83, 12.6, 12.21, 11.73, 11.21, 10.45, 9.49, 8.49, 8.13, 8.19],
                   [6.81, 7.15, 7.68, 8.25, 8.34, 8.09, 7.89, 7.96, 8.02, 8.03, 8.55, 9.35, 10.25, 11.34, 12.03, 12.14, 12.04, 11.72, 11.19, 10.4, 9.5, 8.25, 7.18, 6.86, 6.81],
                   [5.53, 5.6, 5.94, 5.96, 5.73, 5.45, 5.41, 5.33, 5.45, 5.58, 5.93, 6.53, 8.68, 8.75, 10.03, 11.22, 11.31, 11.1, 10.47, 9.62, 8.42, 7.08, 6.22, 5.61, 5.53],
                   [4.42, 4.27, 4.34, 4.39, 4.17, 3.97, 3.6, 3.52, 3.56, 3.74, 4.22, 4.84, 6.02, 7.21, 8.43, 9.01, 10.42, 10.27, 9.66, 8.7, 7.3, 6.17, 5.2, 4.49, 4.42],
                   [3.55, 3.42, 3.41, 3.38, 2.9, 2.57, 2.26, 2.12, 2.19, 2.3, 2.69, 3.29, 4.53, 5.17, 6.06, 7.75, 9.07, 9.23, 8.68, 7.62, 6.57, 5.46, 4.34, 3.68, 3.55],
                   [2.87, 2.78, 2.47, 2.24, 1.96, 1.6, 1.28, 1.22, 1.17, 1.29, 1.55, 2.03, 3.06, 3.94, 4.57, 5.53, 7.23, 8.03, 7.82, 7.02, 5.8, 4.42, 3.56, 3.11, 2.87],
                   [2.32, 2.07, 1.84, 1.59, 1.24, 0.94, 0.71, 0.58, 0.57, 0.66, 0.85, 1.18, 2.02, 2.66, 3.34, 4.17, 5.04, 5.92, 6.28, 5.52, 4.5, 3.6, 2.93, 2.49, 2.32],
                   [1.85, 1.57, 1.36, 1.07, 0.81, 0.54, 0.34, 0.2, 0.22, 0.26, 0.36, 0.61, 1.18, 1.68, 2.3, 3.04, 3.87, 4.37, 4.51, 4.24, 3.64, 2.9, 2.31, 1.92, 1.85],
                   [1.43, 1.17, 0.95, 0.69, 0.49, 0.28, 0.11, 0.05, 0.03, 0.04, 0.12, 0.25, 0.67, 1.04, 1.53, 2.02, 2.62, 3.15, 3.46, 3.24, 2.78, 2.29, 1.94, 1.47, 1.43],
                   [1.06, 0.84, 0.66, 0.46, 0.2, 0.12, 0.03, 0.0, 0.0, 0.0, 0.02, 0.06, 0.34, 0.59, 0.92, 1.32, 1.73, 2.06, 2.34, 2.31, 2.02, 1.7, 1.39, 1.14, 1.06],
                   [0.76, 0.55, 0.39, 0.24, 0.09, 0.03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.11, 0.32, 0.51, 0.78, 1.03, 1.24, 1.42, 1.46, 1.37, 1.16, 1.0, 0.79, 0.76],
                   [0.48, 0.35, 0.22, 0.09, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03, 0.17, 0.24, 0.43, 0.6, 0.72, 0.8, 0.84, 0.85, 0.72, 0.65, 0.51, 0.48],
                   [0.26, 0.17, 0.12, 0.06, 0.03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03, 0.05, 0.09, 0.2, 0.32, 0.38, 0.44, 0.43, 0.44, 0.41, 0.37, 0.29, 0.26],
                   [0.12, 0.07, 0.06, 0.02, 0.02, 0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03, 0.03, 0.08, 0.07, 0.11, 0.11, 0.18, 0.16, 0.16, 0.16, 0.16, 0.12],
                   [0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]]

neutron_detectors = [["AATA", "Alma-Ata A", 5.9, 897], ["AATB", "Alma-Ata B", 5.9, 3340], ["AHMD", "Ahmedabad", 15.94, 50], ["APTY", "Apatity", 0.65, 181],
                     ["ARNM", "Aragats", 7.1, 3200], ["ATHN", "Athens", 8.53, 260], ["BKSN", "Baksan", 5.7, 1700], ["CALG", "Calgary", 1.08, 1123],
                     ["CALM", "NM de Castilla la Mancha", 6.95, 708], ["CLMX", "Climax", 3.0, 3400], ["DJON", "Daejeon", 11.2, 200],
                     ["DOMB", "Dome C mini NM (bare)", 0.01, 3233], ["DOMC", "Dome C mini NM", 0.01, 3233], ["DRBS", "Dourbes", 3.18, 225],
                     ["ESOI", "Emilio Segre Obs. Israel", 10.75, 2055], ["FSMT", "Forth Smith", 0.3, 180], ["HRMS", "Hermanus", 4.58, 26],
                     ["HUAN", "Huancayo", 12.92, 3400], ["INVK", "Inuvik", 0.3, 21], ["IRK2", "Irkustk 2", 3.64, 2000], ["IRK3", "Irkustk 3", 3.64, 3000],
                     ["IRKT", "Irkustk", 3.64, 435], ["JBGO", "JangBogo", 0.3, 29], ["JUNG", "IGY Jungfraujoch", 4.49, 3570],
                     ["JUNG1", "NM64 Jungfraujoch", 4.49, 3475], ["KERG", "Kerguelen", 1.14, 33], ["KGSN", "Kingston", 1.88, 65],
                     ["KIEL", "Kiel", 2.36, 54], ["KIEL2", "KielRT", 2.36, 54], ["LMKS", "Lomnicky Stit", 3.84, 2634], ["MCMU", "Mc Murdo", 0.3, 48],
                     ["MCRL", "Mobile Cosmic Ray Laboratory", 2.46, 200], ["MGDN", "Magadan", 2.1, 220], ["MOSC", "Moscow", 2.43, 200], ["MRNY", "Mirny", 0.03, 30],
                     ["MWSN", "Mawson", 0.22, 30], ["MXCO", "Mexico", 8.28, 2274], ["NAIN", "Nain", 0.3, 46], ["NANM", "Nor-Amberd", 7.1, 2000],
                     ["NEU3", "Neumayer III mini neutron monitor", 0.1, 40], ["NEWK", "Newark", 2.4, 50], ["NRLK", "Norilsk", 0.63, 0],
                     ["NVBK", "Novosibirsk", 2.91, 163], ["OULU", "Oulu", 0.81, 15], ["PSNM", "Doi Inthanon (Princess Sirindhorn NM)", 16.8, 2565],
                     ["PTFM", "Potchefstroom", 6.98, 1351], ["PWNK", "Peawanuck", 0.3, 53], ["ROME", "Rome", 6.27, 0], ["SANB", "Sanae D", 0.73, 52],
                     ["SNAE", "Sanae IV", 0.73, 856], ["SOPB", "South Pole Bare", 0.1, 2820], ["SOPO", "South Pole", 0.1, 2820], ["TERA", "Terre Adelie", 0.01, 32],
                     ["THUL", "Thule", 0.3, 26], ["TSMB", "Tsumeb", 9.15, 1240], ["TXBY", "Tixie Bay", 0.48, 0], ["UFSZ", "Zugspitze", 4.1, 2650],
                     ["YKTK", "Yakutsk", 1.65, 105], ["ZUGS", "Zugspitze", 4.24, 2960]]
