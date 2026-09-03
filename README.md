# Eisabwurf


#### **Eisabfall-Simulation von Windenergieanlagen**



Dieses Python-Programm wurde im Rahmen einer Bachelorarbeit zur numerischen Untersuchung der Flugweite

eines idealisierten Eisstücks beim Eisabfall von Windenergieanlagen entwickelt.



Die Anwendung berechnet die zweidimensionale Flugbahn eines Eisobjekts unter Berücksichtigung der Gravitation,

des aerodynamischen Widerstands, der Windgeschwindigkeit sowie der tangentialen Startgeschwindigkeit.



Zur numerischen Berechnung können folgende Verfahren verwendet und miteinander Verglichen werden:

* Explizites Eulerverfahren
* Runge-Kutta-4-Verfahren
* Verlet-Verfahren



Die Anwendung verfügt über eine grafische Benutzeroberfläche auf Basis von Streamlit. Sämtliche 

relevante Modellparameter können über die Sidebar angepasst werden. Die im Programm hinterlegten Standartwerte

entsprechen den Parametern, die in der BA angesetzt wurden. 



Als Ergebnis werden insbesondere die berechneten Flugweiten der ausgewählten Verfahren, die Flugbahnen des 

Eisobjekts sowie der maximale simulierte Gefährdungsradius ausgegeben. Der Gefährdungsradius wird

auf einer interaktiven Karte dargestellt (am beispielhaft gewählten Standort Cuxhaven).



**Voraussetzungen**



Zur Ausführung werden Python 3 sowie folgende Python-Bibliotheken benötigt:

* streamlit
* folium
* streamlit-folium
* numpy
* pandas



Die benötigten Bibliotheken können beispielsweise über pip installiert werden:

pip install streamlit folium streamlit-folium numpy Pandas



**Ausführung**



Die Anwendung wird über ein Terminal im Verzeichnis der Python-Datei gestartet:

python -m streamlit run eisabfall\_simulation.py (siehe letzte Code-Zeile)



**Einstellbare Parameter**



Über die Benutzeroberfälche können folgende Parameter verändert werden:

* geografische Position der Windenergieanlage
* Nabenhöhe und Rotorradius
* Referenzwindgeschwindigkeit und Messhöhe
* Exponent des Hellmann-Windprofils
* Masse des Eisobjekts
* projizierte Fläche des Eisobjekts
* Luftwiderstandsbeiwert
* Trudeldrehzahl der Windenergieanlage
* numerischer Zeitschritt
* Windprofil-Modellierung
* verwendetes numerisches Lösungsverfahren
