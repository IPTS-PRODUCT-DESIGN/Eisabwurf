import streamlit as st # Benutzeroberfläche
import folium # Kartendarstellung
from streamlit_folium import st_folium # (Die Brücke)
import numpy as np # Numerische Berechnungen
import pandas as pd # Datenverarbeitung mit DataFrames

st.set_page_config(page_title="Eisabfall", layout="wide")
st.title("Eisabfall-Simulation: Numerische Flugbahn")

# =============================================================================
# 1. SIDEBAR PARAMETER 
# =============================================================================
st.sidebar.header("Geografische Position der WEA")
lat0 = st.sidebar.number_input("Breitengrad (Latitude)", value=53.822029, format="%.6f")
lon0 = st.sidebar.number_input("Längengrad (Longitude)", value=8.731000, format="%.6f")

st.sidebar.header("1. Geometrie der WEA")
h_nabe = st.sidebar.number_input("Nabenhöhe H (m)", value=174.50)
r_rotor = st.sidebar.number_input("Rotorradius R (m)", value=87.50)
h_max = h_nabe + r_rotor
st.sidebar.info(f"Maximale Höhe der Blattspitze: {h_max} m")

st.sidebar.header("2. Wind- & Klimadaten")
h_mess = st.sidebar.number_input("Messmasthöhe für Referenzwind (m)", value=21.924)
v_999 = st.sidebar.number_input("99.9% Quantil Windgeschwindigkeit (m/s)", value=17.60)
alpha = st.sidebar.slider("Windprofil-Exponent (alpha)", 0.10, 0.30, 0.16, step=0.01)

st.sidebar.header("3. Aerodynamik des Eisobjekts")
st.sidebar.caption("Modell: Dünne, flache Platte (Konservativer Ansatz)")
m_eis = st.sidebar.number_input("Eisobjekt-Masse (kg)", value=1.0, format="%.3f")
a_proj = st.sidebar.number_input("Projektionsfläche A (m²)", value=0.0347052, format="%.7f")
cw_wert = st.sidebar.number_input("Luftwiderstandsbeiwert (cw)", value=1.31346, format="%.5f")

st.sidebar.header("4. Betriebszustand")
n_trudel = st.sidebar.slider("Trudeldrehzahl (U/min)", 0.0, 5.0, 1.5, step=0.1)

st.sidebar.header("5. Methodik & Evaluation")
dt = st.sidebar.selectbox("Numerischer Zeitschritt dt (s)", options=[0.01, 0.1, 0.5], index=0)
wind_mode = st.sidebar.radio(
    "Windfeld-Modellierung:",
    options=["Dynamisches Hellmann-Profil", "Konstantes Windfeld (TÜV-Ansatz)"]
)

algos = st.sidebar.multiselect(
    "Zu verwendende Algorithmen:",
    options=["Expliziter Euler", "Runge-Kutta 4 (RK4)", "Verlet-Verfahren"],
    default=["Expliziter Euler", "Runge-Kutta 4 (RK4)", "Verlet-Verfahren"]
)

# =============================================================================
# 2. PHYSIKALISCHE GRUNDLAGEN-FUNKTIONEN
# =============================================================================
g = 9.81      
rho = 1.3   

def get_wind_speed(z_pos):
    if wind_mode == "Dynamisches Hellmann-Profil":
        z_safe = max(z_pos, 0.1) # Vermeidung einer Höhe von z = 0 bei der Berechnung des Windprofils
        return v_999 * ((z_safe / h_mess) ** alpha) 
    else:
        return v_999 * ((h_nabe / h_mess) ** alpha)

def get_accelerations(vx_val, vz_val, z_val): 
    v_w_local = get_wind_speed(z_val)
    v_rel_x = v_w_local - vx_val # Relative Geschwindigkeit in x-Richtung
    v_rel_z = 0.0 - vz_val # Relative Geschwindigkeit in z-Richtung 
    v_rel_norm = np.sqrt(v_rel_x**2 + v_rel_z**2) # Normale Geschwindigkeit der Relativströmung (über Satz des Pythagoras)
    
    F_w_factor = 0.5 * rho * cw_wert * a_proj * v_rel_norm # Gemeinsamer Faktor der aerodynamischen Widerstandskraft(klassische Strömungsmechanik)
    
    # Vorzeichenkorrektur: Die Kraft wirkt IMMER in Richtung der Relativströmung / die Beschleunigungen nach Newtons zweitem Gesetz (F = m * a) berechnen
    ax_val = (F_w_factor * v_rel_x) / m_eis
    az_val = -g + ((F_w_factor * v_rel_z) / m_eis)
    return ax_val, az_val

# Startgeschwindigkeit aus Trudeldrehzahl
v_t = ((2 * np.pi * n_trudel) / 60.0) * r_rotor  # Tangentialgeschwindigkeit der Blattspitze in m/s

# =============================================================================
# 3. DIE NUMERISCHEN LOOP-ALGORITHMEN
# =============================================================================

def run_euler():
    x, z = 0.0, h_max
    vx, vz = v_t, 0.0
    path_x, path_z = [x], [z]
    while z > 0:
        ax, az = get_accelerations(vx, vz, z)
        x_next = x + vx * dt
        z_next = z + vz * dt
        
        if z_next <= 0:
            fraction = z / (z - z_next) #lineare Interpolation, um den genauen Aufprallpunkt zu bestimmen (z=0)
            x += vx * dt * fraction
            path_x.append(x)
            path_z.append(0.0)
            break
            
        x, z = x_next, z_next
        vx += ax * dt
        vz += az * dt
        path_x.append(x)
        path_z.append(z)
    return path_x, path_z, x

def run_rk4():
    x, z = 0.0, h_max
    vx, vz = v_t, 0.0
    path_x, path_z = [x], [z]
    while z > 0:
        # k1 (Startpunkt)
        ax1, az1 = get_accelerations(vx, vz, z)
        
        # k2 (Mitte des Intervalls mit k1-Vorschau)
        x2 = x + 0.5 * dt * vx
        z2 = z + 0.5 * dt * vz
        vx2 = vx + 0.5 * dt * ax1
        vz2 = vz + 0.5 * dt * az1
        ax2, az2 = get_accelerations(vx2, vz2, z2)
        
        # k3 (Mitte des Intervalls mit k2-Vorschau)
        x3 = x + 0.5 * dt * vx2
        z3 = z + 0.5 * dt * vz2
        vx3 = vx + 0.5 * dt * ax2
        vz3 = vz + 0.5 * dt * az2
        ax3, az3 = get_accelerations(vx3, vz3, z3)
        
        # k4 (Endpunkt des Intervalls mit k3-Vorschau)
        x4 = x + dt * vx3
        z4 = z + dt * vz3
        vx4 = vx + dt * ax3
        vz4 = vz + dt * az3
        ax4, az4 = get_accelerations(vx4, vz4, z4)
        
        # Gewichtete Updates für Position und Geschwindigkeit
        dx = (dt / 6.0) * (vx + 2*vx2 + 2*vx3 + vx4)
        dz = (dt / 6.0) * (vz + 2*vz2 + 2*vz3 + vz4)
        
        if z + dz <= 0:
            fraction = z / (z - (z + dz))
            x += dx * fraction
            path_x.append(x)
            path_z.append(0.0)
            break
            
        x += dx
        z += dz
        vx += (dt / 6.0) * (ax1 + 2*ax2 + 2*ax3 + ax4)
        vz += (dt / 6.0) * (az1 + 2*az2 + 2*az3 + az4)
        path_x.append(x)
        path_z.append(z)
    return path_x, path_z, x

def run_verlet(): 
    vx_init, vz_init = v_t, 0.0
    ax0, az0 = get_accelerations(vx_init, vz_init, h_max)
    x_prev = 0.0 - vx_init * dt + 0.5 * ax0 * (dt**2) # Berechnung des vorherigen Zeitschritts mittels Taylor-Entwicklung
    z_prev = h_max - vz_init * dt + 0.5 * az0 * (dt**2)
    
    x_curr, z_curr = 0.0, h_max
    path_x, path_z = [x_curr], [z_curr]
    
    while z_curr > 0:
        v_x_est = (x_curr - x_prev) / (2 * dt) # Näherung der Geschwindigkeit aus den Positionen des aktuellen und vorherigen Zeitschritts
        v_z_est = (z_curr - z_prev) / (2 * dt)
        
        ax, az = get_accelerations(v_x_est, v_z_est, z_curr)
        
        x_next = 2 * x_curr - x_prev + ax * (dt**2) # Die neue Position wird dann rein aus den Orten und der Beschleunigung berechnet, die Geschwindigkeit wird nicht direkt integriert.
        z_next = 2 * z_curr - z_prev + az * (dt**2)
        
        if z_next <= 0:
            fraction = z_curr / (z_curr - z_next)
            x_curr += (x_next - x_curr) * fraction
            path_x.append(x_curr)
            path_z.append(0.0)
            break
            
        x_prev, z_prev = x_curr, z_curr
        x_curr, z_curr = x_next, z_next
        path_x.append(x_curr)
        path_z.append(z_curr)
    return path_x, path_z, x_curr

# =============================================================================
# 4. EVALUATION & DATA-FRAMING
# =============================================================================
max_wurfweite = 0.0
metrics = {}
plot_dict = {}

if "Runge-Kutta 4 (RK4)" in algos:
    px, pz, mx = run_rk4()
    metrics["Runge-Kutta 4 (RK4) Wurfweite"] = mx
    max_wurfweite = max(max_wurfweite, mx)
    plot_dict["Runge-Kutta 4 (RK4)"] = pd.Series(pz, index=px) #Die x- und z-Koordinaten werden in eine Pandas-Serie gepackt.

if "Expliziter Euler" in algos:
    px, pz, mx = run_euler()
    metrics["Expliziter Euler Wurfweite"] = mx
    max_wurfweite = max(max_wurfweite, mx)
    plot_dict["Expliziter Euler"] = pd.Series(pz, index=px)

if "Verlet-Verfahren" in algos:
    px, pz, mx = run_verlet()
    metrics["Verlet-Verfahren Wurfweite"] = mx
    max_wurfweite = max(max_wurfweite, mx)
    plot_dict["Verlet-Verfahren"] = pd.Series(pz, index=px)

if max_wurfweite == 0.0:
    max_wurfweite = 10.0 # Sicherheitswert für die Kartendarstellung, falls kein Algorithmus ausgewählt wurde

if plot_dict:
    df_chart = pd.DataFrame(plot_dict).sort_index().interpolate(method="linear") #pd.DataFrame(plot_dict) nimmt die einzelnen Serien und wirft sie in eine große gemeinsame Tabelle...sort_index() sortiert die X-Achse chronologisch von 0 bis zum Einschlag. .interpolate(method="linear") zieht dann mathematisch saubere, gerade Linien zwischen den Lücken.
    
    # Kurven strikt nach ihrem Aufprallpunkt abschneiden
    if "Runge-Kutta 4 (RK4)" in df_chart.columns and "Runge-Kutta 4 (RK4) Wurfweite" in metrics:
        df_chart.loc[df_chart.index > metrics["Runge-Kutta 4 (RK4) Wurfweite"], "Runge-Kutta 4 (RK4)"] = None
    if "Expliziter Euler" in df_chart.columns and "Expliziter Euler Wurfweite" in metrics:
        df_chart.loc[df_chart.index > metrics["Expliziter Euler Wurfweite"], "Expliziter Euler"] = None
    if "Verlet-Verfahren" in df_chart.columns and "Verlet-Verfahren Wurfweite" in metrics:
        df_chart.loc[df_chart.index > metrics["Verlet-Verfahren Wurfweite"], "Verlet-Verfahren"] = None

# =============================================================================
# 5. OBERFLÄCHE & VISUALISIERUNG
# =============================================================================
col1, col2 = st.columns([1, 2])

with col1:
    st.subheader("Physikalische Kennzahlen")
    st.metric("Maximaler simulierter Radius", f"{max_wurfweite:.1f} m")
    for key, value in metrics.items():
        st.metric(key, f"{value:.2f} m")
    st.metric("Startgeschwindigkeit (Tangential)", f"{v_t:.2f} m/s ({v_t*3.6:.1f} km/h)")

with col2:
    st.subheader("Numerische Flugbahn des Eisstücks")
    if algos:
        st.line_chart(df_chart)
    else:
        st.warning("Bitte wähle mindestens einen Algorithmus in der Sidebar aus.")

st.markdown("---") 
st.subheader("Maximaler Gefährdungsradius")

m = folium.Map(location=[lat0, lon0], zoom_start=15.45)
folium.Marker([lat0, lon0], popup="WEA").add_to(m)

folium.Circle(
    location=[lat0, lon0],
    radius=max_wurfweite,
    color='red',
    weight=1.7,
    dash_array='9, 6',
    fill=False,
    popup=f"Numerisch berechnetes Max: {max_wurfweite:.1f}m"
).add_to(m)

# --- FADENKREUZ ---
earth_radius = 6371000.0
cross_length = max(max_wurfweite * 1.2, 300.0)

d_lat = (cross_length / earth_radius) * (180.0 / np.pi) # Ist linear. Die metrische Strecke wird durch den Erdradius geteilt und das resultierende Bogenmaß in Grad umgewandelt (180.0 / np.pi).
d_lon = (cross_length / (earth_radius * np.cos(np.pi * lat0 / 180.0))) * (180.0 / np.pi) # Da die Längenkreise der Erde zu den Polen hin enger werden, musst man durch den Kosinus des aktuellen Breitengrads teilen, um die Erdkrümmung zu kompensieren.

folium.PolyLine(locations=[[lat0, lon0 - d_lon], [lat0, lon0 + d_lon]], color='grey', weight=1.5, opacity=0.8, dash_array='5, 5').add_to(m)
folium.PolyLine(locations=[[lat0 - d_lat, lon0], [lat0 + d_lat, lon0]], color='grey', weight=1.5, opacity=0.8, dash_array='5, 5').add_to(m)

step_meters = 100.0
for dist in np.arange(step_meters, cross_length, step_meters):
    m_lat = (dist / earth_radius) * (180.0 / np.pi)
    m_lon = (dist / (earth_radius * np.cos(np.pi * lat0 / 180.0))) * (180.0 / np.pi)
    
    folium.map.Marker([lat0, lon0 + m_lon], icon=folium.DivIcon(html=f'<div style="font-size: 9px; color: grey; font-weight: bold;">┼ {int(dist)}m</div>')).add_to(m)
    folium.map.Marker([lat0, lon0 - m_lon], icon=folium.DivIcon(html=f'<div style="font-size: 9px; color: grey; font-weight: bold;">┼ {int(dist)}m</div>')).add_to(m)
    folium.map.Marker([lat0 + m_lat, lon0], icon=folium.DivIcon(html=f'<div style="font-size: 9px; color: grey; font-weight: bold;">─ {int(dist)}m</div>')).add_to(m)
    folium.map.Marker([lat0 - m_lat, lon0], icon=folium.DivIcon(html=f'<div style="font-size: 9px; color: grey; font-weight: bold;">─ {int(dist)}m</div>')).add_to(m)

st_folium(m, width=1300, height=800, returned_objects=[])

# Start der Streamlit-Anwendung (VS Code):
# python -m streamlit run c:\BA_Eisabfall\eisabfall_simulation.py