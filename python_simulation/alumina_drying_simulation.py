import kaleido
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.integrate import solve_ivp

def positive_arithmetic_remainder(n, m):
    return n - abs(m)*np.floor(n/abs(m))

def upward_drying(t, Y):
    M =  []
    W =  []
    Tm = []
    dWdy = []
    dTmdy = []
    Pab = []
    Psat = []
    UR = []
    Meq_40 = []
    Meq_50 = []
    Meq_60 =  []
    Meq = []
    Rho = []
    Ep = []
    K = []
    dMdt = []
    f = []

    dY = np.zeros(3*N)
    for i in range(N):
        M.append(Y[i])
        W.append(Y[N+i])
        Tm.append(Y[2*N+i])

        if i == 0:
            dWdy.append((W[i] - W0) / hy) 
            dTmdy.append((Tm[i] - Tf0ac) / hy)
        else:
            dWdy.append((W[i] - W[i-1]) / hy)
            dTmdy.append((Tm[i] - Tm[i-1]) / hy) 

        Pab.append((28.97 / 18*W[i]) / (1+28.97 / 18*W[i]) * (695.1 / 760))
        Psat.append(np.exp(18.3036 - 3816.44 / (Tm[i] + 227.02)) / 760)
        UR.append(Pab[i] / Psat[i])
        # Equilibrium isotherm for alumina at 50°C
        c_iso = 1.39796669; 
        k_iso = 0.14982732; 
        mo_iso = 1.148646602; 
        Meq_60.append(0.0358 * UR[i] / (1 - 0.835 * UR[i])**2)
        Meq_40.append(mo_iso * c_iso * k_iso * UR[i] / ((1 - k_iso * UR[i])*(1 - k_iso * UR[i] + k_iso * c_iso * UR[i])))
        Meq_50.append((Meq_40[i] + Meq_60[i]) / 2) 
        Meq.append(Meq_60[i])
        Rho.append(Rhoo)
        Ep.append(Epsilon)
        Ko = 1.124342559
        K.append(Ko * np.exp(-2129.52871 / (Tm[i] + 273.15)))
        dMdt.append(-K[i] * (M[i] - Meq[i])) 
        f.append(-(1 - Ep[i]) * Rho[i] * dMdt[i]) 
        dY[i] = dMdt[i] # Mass balance for the solid phase
        dY[i+N] = (f[i] - (Gf * dWdy[i]) / (Ep[i] * Rhof)) # Mass balance for the fluid phase
        dY[i+2*N] = (-f[i] * (_lambda) - (Gf * (Cpf + W[i] * Cpv) * dTmdy[i])) / (((1 - Ep[i]) * Rho[i] * (Cps + M[i] * Cpl)) + (Ep[i] * Rhof * (Cpf + W[i]))) # Energy balance for the mixture
    return dY

def downard_drying(t, Y):
    M=W=Tm=dWdy=dTmdy=Pab=Psat=UR=Meq_40=Meq_50=Meq_60=Meq=Rho=Ep=K=dMdt=f=np.zeros(N)

    dY = np.zeros(3*N)
    for i in range(N-1,-1,-1):
        M = Y[i]
        W = Y[N+i]
        Tm = Y[2*N+i]

        if i == N-1:
            dWdy = (W - W0) / hy
            dTmdy = (Tm - Tf0des) / hy
        else:
            dWdy = (W - last_W) / hy
            dTmdy = (Tm - last_Tm) / hy

        Pab = (28.97 / 18*W) / (1+28.97 / 18*W) * (695.1 / 760)
        Psat = (np.exp(18.3036 - 3816.44 / (Tm + 227.02)) / 760)
        UR = Pab / Psat
        # Equilibrium isotherm for alumina at 50°C 
        c_iso = 1.39796669
        k_iso = 0.14982732 
        mo_iso = 1.148646602 
        Meq_60 = 0.0358 * UR / (1 - 0.835 * UR)**2
        Meq_40 = mo_iso * c_iso * k_iso * UR / ((1 - k_iso * UR)*(1 - k_iso * UR + k_iso * c_iso * UR))
        Meq_50 = (Meq_40 + Meq_60) / 2
        Meq = Meq_60
        Rho = Rhoo 
        Ep = Epsilon 
        Ko = 1.124342559
        K = Ko * np.exp(-2129.52871 / (Tm + 273.15))
        dMdt = -K * (M - Meq)
        f = -(1 - Ep) * Rho * dMdt
        dY[i] = dMdt # Mass balance for the solid phase
        dY[i+N] = (f - (Gf * dWdy) / (Ep * Rhof)) # Mass balance for the fluid phase
        dY[i+2*N] = (-f * (_lambda) - (Gf * (Cpf + W * Cpv) * dTmdy)) / (((1 - Ep) * Rho * (Cps + M * Cpl)) + (Ep * Rhof * (Cpf + W))) # Energy balance for the mixture

        last_W = W
        last_Tm = Tm
    return dY

def create_df_from_simulation_results(df):
    return pd.DataFrame(
        data = df,
        columns = [
            'Tempo_min',
            'X_0cm', 
            'X_1cm', 
            'X_2cm', 
            'X_3cm', 
            'X_4cm', 
            'X_5cm', 
            'X_6cm', 
            'X_7cm', 
            'X_8cm', 
            'X_9cm', 
            'X_10cm', 
            'H_0cm', 
            'H_1cm',
            'H_2cm',
            'H_3cm',
            'H_4cm',
            'H_5cm',
            'H_6cm',
            'H_7cm',
            'H_8cm',
            'H_9cm',
            'H_10cm',
            'T_0cm',
            'T_1cm',
            'T_2cm',
            'T_3cm',
            'T_4cm',
            'T_5cm',
            'T_6cm',
            'T_7cm',
            'T_8cm',
            'T_9cm',
            'T_10cm'
        ]
    )

def plot_or_save_results(df, plot_type, filename, action = 'plot'):
    plot_title = 'Umidade da alumina em função do tempo'
    yaxis_title = 'Umidade da partícula em base seca [kg água/kg sólido seco]'
    legend_x_position = .83
    legend_y_position = 1
    if plot_type == 'T':
        plot_title = 'Temperatura da mistura em função do tempo'
        yaxis_title = 'Temperatura da mistura [°C]'
        legend_x_position = .81
        legend_y_position = .18    

    fig = go.Figure(
        data = [
            go.Scatter(
                x = df['Tempo_min'],
                y = df[f'{plot_type}_1cm'],
                name = f'{plot_type}(1cm)',
                marker = {
                    'color': plots_series_colors['1'],
                    'line_color': 'black',
                    'line_width': 1.5,
                    
                },
            ),
            go.Scatter(
                x = df['Tempo_min'],
                y = df[f'{plot_type}_5cm'],
                name = f'{plot_type}(5cm)',
                marker = {
                    'color': plots_series_colors['5'],
                    'line_color': 'black',
                    'line_width': 1.5,
                    
                },
            ),
            go.Scatter(
                x = df['Tempo_min'],
                y = df[f'{plot_type}_10cm'],
                name = f'{plot_type}(10cm)',
                marker = {
                    'color': plots_series_colors['10'],
                    'line_color': 'black',
                    'line_width': 1.5,
                    
                },
            ),
        ],
        layout = go.Layout(
            width = 700,
            height = 520,
            title = {
                'font': {
                    'size': 20,
                },
                'text': plot_title,
                'x': .5,
                'y': .87,
            },
            xaxis = {
            'title': {
                'text': 'Tempo [minutos]',
                'font_size': 14,
                }
            },
            yaxis = {
                'title': {
                    'text': yaxis_title,
                    'font_size': 14
                }
            },
            font = {
                'family': 'Times New Roman',
                'color': 'black'
            },
            legend = {
                'x': legend_x_position,
                'y': legend_y_position,
                'bgcolor': '#e5e5e5',
            },
            plot_bgcolor = '#e5e5e5'
        )
    )

    fig.update_xaxes(
        tickmode = 'linear',
        tick0 = 0, 
        dtick = 10
    )
    fig.update_yaxes(
        tickmode = 'linear',
        tick0 = 0,
        dtick = .05
    )
    if plot_type == 'T':
        fig.update_yaxes(
            tickmode = 'linear',
            tick0 = 20,
            dtick = 5
        )

    if action == 'plot':
        fig.show()
    elif action == 'save':
        fig.write_image(filename, format='png', engine='kaleido')


# Plots series hex colors
plots_series_colors = {
    '1': '#bc4b51',
    '5': '#5b8e7d',
    '10': '#f4a259',
}

# Model parameters 
Epsilon = 0.4 # Porosity of the bed
Rhoo = 1.69 # Alumina specific mass (g/cm**3) 
Cpf = 0.25 # Dry air specific heat (cal/g.°C) 
Cpv = 0.28 # Water vapor air specific heat (cal/g.°C) 
Cps = 0.199914 # Solid specific heat (cal/g.°C) 
Cpl = 1.0 # Liquid water specific heat (cal/g.°C) 
_lambda = 573 # Latent heat of water vaporization (cal/g) 
bed_height = 10.0 # Bed height (cm) 
m = 0.390 # Mass flow rate of air (kg/min) 
D = 10.0 # Bed diameter (cm) 

# Air reversal parameters
t0rev = 10*60 # Time when the first airflow reversal happens
deltrev = 10*60 # Interval between airflow reversals
nrev = 11 # Number of airflow reversals

# Numerical method parameters
N = 11 # Number of reversions
hyaux = np.linspace(0, bed_height, N)
hy = hyaux[1] - hyaux[0]

# Experiment conditions
Ufo = 0.016 # Initial air humidity
Uso = 0.45 # Initial alumina humidity
Tfo = 60.0 # Initial bed temperature
Tf0ac = Tfo # Initial inlet temperature of the upward flow [°C] 
Tf0des = Tfo # Initial inlet temperature of the downward flow [°C] 
Tmo = 20.8 # Initial bed temperature, mixture of fluid and solid phases [°C] 

W0 = Ufo
M0 = Uso 
Rhof = (0.0012*293.15 / (Tfo+273.15)) * (1 / (1+Ufo)) # Fluid density (g/cm**3) 
Gf = m*1000 / (60*np.pi*(D**2) / 4) # Mass flux of the fluid phase [g/cm**2.s]
t0 = 0.001 # Initial drying time (s)
delt = 20 # Simulation time step (s)
tf = t0rev # Final drying time (s)

y0 = np.array([Uso*np.ones(N), Ufo*np.ones(N), Tmo*np.ones(N)]).flatten() # Simulation initial state
simulation_results = y0.reshape((1, len(y0))) # Array that stores the results of the simulation
time_steps = [0] # List that stores the time steps of the simulation

for i in range(N):
    if positive_arithmetic_remainder(i,2) == 0: # Upward flow
        t_eval = np.arange(t0, tf, delt)
        t_span = [t0, tf]
        if i > 0:
            t_eval = np.arange(t0, tf+.05, delt)
            t_span = [t0, tf+.05]
        Y = solve_ivp(
            fun = upward_drying,
            t_span = t_span,
            t_eval = t_eval,
            y0 = y0,
            method = 'Radau',
            atol = 1e-10,
            rtol = 1e-10
        )
        y0 = np.transpose(Y.y)[-1]
        t0 = tf
        tf = tf + deltrev
    elif positive_arithmetic_remainder(i,2) == 1: # Downward flow
        t_eval = np.arange(t0, tf+.05, delt) 
        Y = solve_ivp(
            fun = downard_drying,
            t_span = [t0, tf+.05],
            t_eval = t_eval,
            y0 = y0,
            method = 'Radau',
            atol = 1e-10,
            rtol = 1e-10
        )
        y0 = Y.y[:,-1]
        t0 = tf
        tf = tf + deltrev
    
    simulation_results = np.append(simulation_results, np.transpose(Y.y), axis = 0)
    time_steps.extend(t_eval)

time_steps = [t/60 for t in time_steps] # Time steps values in minutes
time_steps_simulation_results = np.hstack((np.asarray(time_steps).reshape(len(time_steps),1),simulation_results))
time_steps_simulation_results_df = create_df_from_simulation_results(time_steps_simulation_results)

plot_or_save_results(
    df = time_steps_simulation_results_df,
    plot_type = 'X',
    action = 'plot',
    filename = 'simulations_results/python_simulation_results.png'
)

time_steps_simulation_results_df.to_csv(
    'simulations_results/simulation_results_python.csv',
    index=False
)
