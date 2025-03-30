import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

# Enhanced apheresis system parameters with hematocrit factors
APHERESIS_SETTINGS = {
    'Spectra Optia': {
        'interface_range': (0.5, 2.0),
        'flow_range': (40, 70),
        'plasma_removal_range': (5, 25),
        'acd_ratio_range': (12, 16),
        'hct_impact': 0.3,  # Hematocrit impact factor (lower for Optia)
        'rbc_sensitivity': 0.4  # RBC contamination sensitivity
    },
    'Haemonetics': {
        'interface_range': (0.5, 2.0),
        'flow_range': (40, 65),
        'plasma_removal_range': (5, 20),
        'acd_ratio_range': (12, 15),
        'hct_impact': 0.7,  # Higher hematocrit impact
        'rbc_sensitivity': 0.8  # More sensitive to RBCs
    }
}

# UV-C bag parameters (254nm)
BAG_TYPES = {
    'Spectra Optia (Polyethylene)': {'absorption': 0.9, 'scattering': 5.5, 'thickness': 0.20},
    'Haemonetics (PVC)': {'absorption': 1.3, 'scattering': 7.0, 'thickness': 0.25}
}

def calculate_ecp_uvc(tlc, lymph_percent, hct, system, lamp_power, target_dose, 
                     use_hood, bag_type, interface_pos, flow_rate, 
                     plasma_removal, acd_ratio):
    """Enhanced ECP UV-C calculator with hematocrit-adjusted apheresis parameters"""
    
    params = APHERESIS_SETTINGS[system]
    
    # 1. Hematocrit-adjusted apheresis performance
    hct_factor = 1 - (params['hct_impact'] * (hct - 40)/40  # Normalized to 40% Hct
    interface_factor = (1 - (interface_pos - params['interface_range'][0]) / 
                      (params['interface_range'][1] - params['interface_range'][0])) * hct_factor
    
    flow_factor = (flow_rate / params['flow_range'][1]) * hct_factor
    purity_factor = 0.7 + (interface_pos * 0.15 * hct_factor)
    
    # 2. Product composition with Hct-adjusted RBC contamination
    mnc_conc = (tlc * (lymph_percent/100) * 1.2 * 4 * flow_factor * interface_factor)
    rbc_contam = np.mean([2.0, 5.0]) * (1 - plasma_removal/20) * (hct/40) * params['rbc_sensitivity']
    
    # 3. UV-C delivery (unchanged)
    transmission = np.exp(-np.sqrt(3 * BAG_TYPES[bag_type]['absorption'] * 
                          (BAG_TYPES[bag_type]['absorption'] + BAG_TYPES[bag_type]['scattering'])) * 
                         BAG_TYPES[bag_type]['thickness'])
    distance = 20 if use_hood else 15
    effective_intensity = (lamp_power * 1000 * 0.85 * transmission) / (4 * np.pi * distance**2)
    
    # 4. Dose adjustment with enhanced shielding calculation
    shielding = (0.01 * mnc_conc) + (0.03 * rbc_contam * (hct/40))
    effective_dose = target_dose * transmission * max(1 - shielding, 0.3)  # Minimum 30% dose
    exp_time = (effective_dose / (effective_intensity / 1000)) / 60
    
    return {
        'mnc_conc': mnc_conc,
        'rbc_contam': rbc_contam,
        'purity_factor': purity_factor,
        'effective_dose': effective_dose,
        'exp_time': exp_time,
        'transmission': transmission,
        'hct_factor': hct_factor,
        'system_params': params
    }

def main():
    st.set_page_config(page_title="Enhanced ECP Calculator", layout="wide")
    st.title("Extracorporeal Photopheresis (ECP) Calculator with Hematocrit Adjustment")
    
    with st.sidebar:
        st.header("Patient Parameters")
        tlc = st.slider("TLC (×10³/µL)", 1.0, 50.0, 8.0, 0.5)
        lymph_percent = st.slider("Lymphocyte %", 5, 90, 30)
        hct = st.slider("Hematocrit (%)", 20.0, 60.0, 40.0, 0.1,
                       help="Critical for apheresis efficiency calculation")
        
        st.header("System Configuration")
        system = st.selectbox("Apheresis System", list(APHERESIS_SETTINGS.keys()))
        bag_type = st.selectbox("UV-C Bag Type", list(BAG_TYPES.keys()))
        
        st.header("Treatment Parameters")
        lamp_power = st.slider("UV-C Lamp Power (W)", 5, 50, 25)
        target_dose = st.slider("Target Dose (J/cm²)", 0.0, 10.0, 2.5, 0.1)
        use_hood = st.checkbox("Use Laminar Hood", value=True)
        
        st.header("Apheresis Settings")
        interface_pos = st.slider("Interface Position", 
                                APHERESIS_SETTINGS[system]['interface_range'][0], 
                                APHERESIS_SETTINGS[system]['interface_range'][1], 
                                1.0, 0.1)
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
    results = calculate_ecp_uvc(tlc, lymph_percent, hct, system, lamp_power, target_dose,
                              use_hood, bag_type, interface_pos, flow_rate,
                              plasma_removal, acd_ratio)
    
    # Display results
    st.subheader("Procedure Parameters")
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("System Efficiency Factor", f"{results['hct_factor']:.2f}",
                help=f"Hematocrit-adjusted efficiency for {system}")
        st.metric("MNC Concentration", f"{results['mnc_conc']:.1f} ×10⁶/mL")
        st.metric("RBC Contamination", f"{results['rbc_contam']:.1f} ×10⁹",
                help=f"Higher in Haemonetics at Hct >40%: {APHERESIS_SETTINGS[system]['rbc_sensitivity']*100:.0f}% sensitivity")
        
    with col2:
        st.metric("Effective UV Dose", f"{results['effective_dose']:.2f} J/cm²")
        st.metric("Treatment Time", f"{results['exp_time']:.1f} minutes")
        st.metric("T-cell Purity", f"{results['purity_factor']:.2f}")
    
    # System-specific guidance
    st.subheader("System-Specific Optimization")
    if system == 'Haemonetics':
        st.warning("""
        **For High Hematocrit (>45%):**
        - Reduce flow rate by 10-20%
        - Increase plasma removal by 5%
        - Consider lower interface position (0.5-1.0)
        """)
    else:
        st.info("""
        **Spectra Optia Tips:**
        - More tolerant of high Hct (up to 50%)
        - Maintain flow rate 50-60mL/min for optimal yield
        """)
    
    # Create plots
    st.subheader("Treatment Response Curves")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    
    # Dose-response plot
    doses = np.linspace(0, 5, 100)
    ax1.plot(doses, 100*np.exp(-1.0*doses*results['transmission']*results['purity_factor']), 'r-', label='T-cells')
    ax1.plot(doses, 100*np.exp(-0.25*doses*results['transmission']), 'b-', label='Monocytes')
    ax1.axvline(results['effective_dose'], color='k', linestyle='--', label='Selected Dose')
    ax1.set_title(f'Dose-Response (Hct={hct}%)')
    ax1.set_xlabel('UV-C Dose (J/cm²)')
    ax1.set_ylabel('Viability (%)')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Time-response plot
    times = np.linspace(0, max(results['exp_time']*2, 90), 100)
    time_doses = (results['effective_intensity']/1000) * (times * 60)
    ax2.plot(times, 100*np.exp(-1.0*time_doses*results['purity_factor']), 'r-', label='T-cells')
    ax2.plot(times, 100*np.exp(-0.25*time_doses), 'b-', label='Monocytes')
    ax2.axvline(results['exp_time'], color='k', linestyle='--')
    ax2.set_title('Time-Response Curve')
    ax2.set_xlabel('Time (minutes)')
    ax2.legend()
    ax2.grid(alpha=0.3)
    
    st.pyplot(fig)
    
    # Clinical protocol
    st.subheader("Clinical Protocol Recommendations")
    st.markdown(f"""
    **For Hct {hct}% with {system}:**
    1. **Interface:** Maintain {interface_pos:.1f} ± 0.2
    2. **Flow Rate:** {flow_rate} mL/min (adjust ±10% for Hct >45%)
    3. **ACD Ratio:** 1:{acd_ratio} (increase to 1:{min(16, acd_ratio+1)} if RBC >5×10⁹)
    4. **UV Dose:** Target {results['effective_dose']:.2f} J/cm² ± 15%
    
    **Hematocrit-Specific Notes:**
    - Current efficiency factor: {results['hct_factor']:.2f} (1.0 = ideal at 40% Hct)
    - RBC contamination increased by {(hct/40-1)*100:.0f}% vs normal Hct
    """)

if __name__ == "__main__":
    main()
