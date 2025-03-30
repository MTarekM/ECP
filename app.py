import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

# Apheresis system parameters
APHERESIS_SETTINGS = {
    'Spectra Optia': {
        'interface_range': (0.5, 2.0),
        'flow_range': (40, 70),
        'plasma_removal_range': (5, 25),  # %
        'acd_ratio_range': (12, 16)       # 1:12 to 1:16
    },
    'Haemonetics': {
        'interface_range': (0.5, 2.0),
        'flow_range': (40, 65),
        'plasma_removal_range': (5, 20),
        'acd_ratio_range': (12, 15)
    }
}

# UV-C bag parameters (254nm)
BAG_TYPES = {
    'Spectra Optia (Polyethylene)': {'absorption': 0.9, 'scattering': 5.5, 'thickness': 0.20},
    'Haemonetics (PVC)': {'absorption': 1.3, 'scattering': 7.0, 'thickness': 0.25}
}

def calculate_ecp_uvc(tlc, lymph_percent, system, lamp_power, target_dose, 
                     use_hood, bag_type, interface_pos, flow_rate, 
                     plasma_removal, acd_ratio):
    """Enhanced ECP UV-C calculator with apheresis parameters"""
    
    # 1. Calculate apheresis performance factors
    params = APHERESIS_SETTINGS[system]
    interface_factor = 1 - (interface_pos - params['interface_range'][0]) / (params['interface_range'][1] - params['interface_range'][0])
    flow_factor = flow_rate / params['flow_range'][1]
    purity_factor = 0.7 + (interface_pos * 0.15)  # Higher interface → more T-cells
    
    # 2. Estimate product composition
    mnc_conc = (tlc * (lymph_percent/100) * 1.2 * 4 * flow_factor * interface_factor)
    rbc_contam = np.mean([2.0, 5.0]) * (1 - plasma_removal/20)  # Higher plasma removal → less RBC
    
    # 3. UV-C delivery calculations
    transmission = np.exp(-np.sqrt(3 * BAG_TYPES[bag_type]['absorption'] * 
                                  (BAG_TYPES[bag_type]['absorption'] + BAG_TYPES[bag_type]['scattering'])) * 
                         BAG_TYPES[bag_type]['thickness'])
    distance = 20 if use_hood else 15
    effective_intensity = (lamp_power * 1000 * 0.85 * transmission) / (4 * np.pi * distance**2)
    
    # 4. Dose adjustment with apheresis factors
    shielding = (0.01 * mnc_conc) + (0.02 * rbc_contam)
    effective_dose = target_dose * transmission * max(1 - shielding, 0.5)
    exp_time = (effective_dose / (effective_intensity / 1000)) / 60
    
    return {
        'mnc_conc': mnc_conc,
        'rbc_contam': rbc_contam,
        'purity_factor': purity_factor,
        'effective_dose': effective_dose,
        'exp_time': exp_time,
        'transmission': transmission,
        'effective_intensity': effective_intensity
    }

def main():
    st.set_page_config(page_title="ECP UV-C Calculator", layout="wide")
    st.title("Extracorporeal Photopheresis (ECP) UV-C Calculator")
    
    with st.sidebar:
        st.header("Apheresis Parameters")
        system = st.selectbox("System", list(APHERESIS_SETTINGS.keys()))
        tlc = st.slider("TLC (×10³/µL)", 1.0, 50.0, 8.0, 0.5)
        lymph_percent = st.slider("Lymphocyte %", 5, 90, 30)
        
        st.header("UV-C Treatment Parameters")
        lamp_power = st.slider("UV-C Lamp Power (W)", 5, 30, 15)
        target_dose = st.slider("Target Dose (J/cm²)", 0.0, 6.0, 1.0, 0.1)
        use_hood = st.checkbox("Use Laminar Hood", value=True)
        bag_type = st.selectbox("Bag Type", list(BAG_TYPES.keys()))
        
        st.header("Apheresis Settings")
        interface_pos = st.slider("Interface Position", 0.5, 2.0, 1.0, 0.1)
        flow_rate = st.slider("Flow Rate (mL/min)", 
                             APHERESIS_SETTINGS[system]['flow_range'][0], 
                             APHERESIS_SETTINGS[system]['flow_range'][1], 50)
        plasma_removal = st.slider("Plasma Removal (%)", 
                                  APHERESIS_SETTINGS[system]['plasma_removal_range'][0], 
                                  APHERESIS_SETTINGS[system]['plasma_removal_range'][1], 15)
        acd_ratio = st.slider("ACD Ratio (1:X)", 
                             APHERESIS_SETTINGS[system]['acd_ratio_range'][0], 
                             APHERESIS_SETTINGS[system]['acd_ratio_range'][1], 13)
    
    # Calculate results
    results = calculate_ecp_uvc(tlc, lymph_percent, system, lamp_power, target_dose, 
                               use_hood, bag_type, interface_pos, flow_rate, 
                               plasma_removal, acd_ratio)
    
    # Display results
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Product Characteristics")
        st.metric("Estimated MNC Concentration", f"{results['mnc_conc']:.1f} ×10⁶/mL")
        st.metric("RBC Contamination", f"{results['rbc_contam']:.1f} ×10⁹")
        st.metric("T-cell Purity Factor", f"{results['purity_factor']:.2f}")
        
    with col2:
        st.subheader("UV-C Treatment")
        st.metric("Effective Dose", f"{results['effective_dose']:.2f} J/cm²")
        st.metric("Estimated Treatment Time", f"{results['exp_time']:.1f} minutes")
        st.metric("UV Transmission", f"{results['transmission']*100:.1f}%")
    
    # Create plots
    st.subheader("Dose-Response Curves")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    
    # Dose-response plot
    doses = np.linspace(0, 2.5, 100)
    ax1.plot(doses, 100*np.exp(-1.0*doses*results['transmission']*results['purity_factor']), 'r-', label='T-cells')
    ax1.plot(doses, 100*np.exp(-0.25*doses*results['transmission']), 'b-', label='Monocytes')
    ax1.axvline(results['effective_dose'], color='k', linestyle='--')
    ax1.set_title(f'ECP Dose-Response\nInterface={interface_pos}, Flow={flow_rate}mL/min')
    ax1.set_xlabel('UV-C Dose (J/cm²)')
    ax1.set_ylabel('Viability (%)')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Time-response plot
    times = np.linspace(0, max(results['exp_time']*2, 40), 100)
    time_doses = (results['effective_intensity']/1000) * (times * 60)
    ax2.plot(times, 100*np.exp(-1.0*time_doses*results['purity_factor']), 'r-')
    ax2.plot(times, 100*np.exp(-0.25*time_doses), 'b-')
    ax2.axvline(results['exp_time'], color='k', linestyle='--')
    ax2.set_title('Treatment Time Course')
    ax2.set_xlabel('Time (minutes)')
    ax2.grid(alpha=0.3)
    
    st.pyplot(fig)
    
    # Clinical guidance
    st.subheader("Clinical Guidance")
    st.markdown("""
    - **Lower interface position** → More T-cell depletion
    - **Higher flow rate** → Faster processing but less purity
    - **Adjust plasma removal** to control RBC contamination
    - **Monitor treatment time** based on calculated effective dose
    """)

if __name__ == "__main__":
    main()
