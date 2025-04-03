import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

# Enhanced apheresis system parameters with hematocrit factors
APHERESIS_SETTINGS = {
    'Spectra Optia': {
        'flow_range': (40, 70),
        'plasma_removal_range': (5, 25),
        'acd_ratio_range': (12, 16),
        'hct_impact': 0.3,  # 30% hematocrit sensitivity
        'rbc_base_contam': 2.0  # Base RBC contamination (×10⁹)
    },
    'Haemonetics': {
        'flow_range': (40, 65),
        'plasma_removal_range': (5, 20),
        'acd_ratio_range': (12, 15),
        'hct_impact': 0.7,  # 70% hematocrit sensitivity
        'rbc_base_contam': 5.0  # Higher base contamination
    }
}

BAG_TYPES = {
    'Spectra Optia (Polyethylene)': {'absorption': 0.9, 'scattering': 5.5, 'thickness': 0.20},
    'Haemonetics (PVC)': {'absorption': 1.3, 'scattering': 7.0, 'thickness': 0.25}
}

def calculate_ecp_uv(tlc, lymph_percent, hct, system, lamp_power, target_dose, 
                     use_hood, custom_distance, bag_type, flow_rate, plasma_removal, acd_ratio):
    """Hematocrit-adjusted ECP calculator with system-specific efficiencies"""
    
    params = APHERESIS_SETTINGS[system]
    
    # 1. Hematocrit efficiency correction (normalized to 40% Hct)
    hct_efficiency = 1 - params['hct_impact'] * (hct - 40)/40
    
    # 2. Adjusted apheresis performance with optimal interface position (1.25)
    interface_factor = hct_efficiency  # Optimal at fixed 1.25 position
    flow_factor = (flow_rate - params['flow_range'][0]) / (params['flow_range'][1] - params['flow_range'][0]) * hct_efficiency
    purity_factor = 0.7 + (1.25 * 0.15 * hct_efficiency)  # Interface position fixed at 1.25
    
    # 3. Product composition with Hct-adjusted RBC contamination
    mnc_conc = (tlc * (lymph_percent/100) * 4 * flow_factor * interface_factor)
    
    # RBC contamination increases with Hct and differs by system
    rbc_contam = params['rbc_base_contam'] * (hct/40) * (1 - plasma_removal/25)
    
    # 4. UV delivery calculations
    transmission = np.exp(-np.sqrt(3 * BAG_TYPES[bag_type]['absorption'] * 
                          (BAG_TYPES[bag_type]['absorption'] + BAG_TYPES[bag_type]['scattering'])) * 
                         BAG_TYPES[bag_type]['thickness'])
    distance = custom_distance if not use_hood else 20  # Hood always uses 20cm distance
    effective_intensity = (lamp_power * 1000 * 0.85 * transmission) / (4 * np.pi * distance**2)
    
    # 5. Dose adjustment with Hct-impacted shielding
    shielding = (0.01 * mnc_conc) + (0.05 * rbc_contam * (hct/40))
    effective_dose = target_dose * transmission * max(1 - shielding, 0.3)
    exp_time = (effective_dose / (effective_intensity / 1000)) / 60
    
    return {
        'mnc_conc': mnc_conc,
        'rbc_contam': rbc_contam,
        'purity_factor': purity_factor,
        'effective_dose': effective_dose,
        'exp_time': exp_time,
        'transmission': transmission,
        'effective_intensity': effective_intensity,
        'hct_efficiency': hct_efficiency,
        'system': system,
        'distance': distance
    }

def main():
    st.set_page_config(page_title="UV-based Sensitizer-free ECP Calculator", layout="wide")
    st.title("UV-based Sensitizer-free ECP Calculator")
    
    with st.sidebar:
        st.header("Patient Blood Parameters")
        col1, col2 = st.columns(2)
        with col1:
            tlc = st.number_input("TLC (×10³/µL)", min_value=1.0, max_value=50.0, value=8.0, step=0.5)
        with col2:
            lymph_percent = st.number_input("Lymphocyte %", min_value=5, max_value=90, value=30)
        
        hct = st.slider("Patient Hematocrit (%)", min_value=20.0, max_value=60.0, value=40.0, step=0.1,
                       help="Critical for apheresis efficiency and RBC contamination")
        
        st.header("Treatment Setup")
        system = st.selectbox("Apheresis System", list(APHERESIS_SETTINGS.keys()))
        bag_type = st.selectbox("UV Bag Type", list(BAG_TYPES.keys()))
        
        st.header("UV Parameters")
        uv_type = st.selectbox("UV Type", ["UV-A", "UV-B", "UV-C"], index=2)
        
        if uv_type == "UV-A":
            target_dose = st.slider("Target Dose (J/cm²)", 5.0, 10.0, 5.0, 0.1)
        elif uv_type == "UV-B":
            target_dose = st.slider("Target Dose (J/cm²)", 0.5, 2.0, 1.0, 0.1)
        else:
            target_dose = st.slider("Target Dose (J/cm²)", 2.0, 6.0, 2.5, 0.1)
        
        lamp_power = st.slider("Lamp Power (W)", 5, 50, 25)
        use_hood = st.checkbox("Use Laminar Hood (fixed 20cm distance)", value=True)
        
        if not use_hood:
            custom_distance = st.slider("Custom Distance (cm)", 10, 50, 15, 1,
                                       help="Distance between UV source and treatment bag")
        else:
            custom_distance = 20  # Default when hood is used

    # Apheresis Settings Section with HCT Warnings
    with st.sidebar:
        st.header("Apheresis Settings")
        
        # Flow rate adjustment guidance
        flow_default = 50
        if hct > 45:
            flow_default = 45 if system == 'Haemonetics' else 55
        flow_rate = st.slider("Flow Rate (mL/min)", 
                            APHERESIS_SETTINGS[system]['flow_range'][0], 
                            APHERESIS_SETTINGS[system]['flow_range'][1], 
                            flow_default)
        
        plasma_removal = st.slider("Plasma Removal (%)", 
                                 APHERESIS_SETTINGS[system]['plasma_removal_range'][0], 
                                 APHERESIS_SETTINGS[system]['plasma_removal_range'][1], 15)
        
        # ACD ratio adjustment for high Hct
        acd_default = 13
        if hct > 45:
            acd_default = 12 if system == 'Haemonetics' else 14
        acd_ratio = st.slider("ACD Ratio (1:X)", 
                            APHERESIS_SETTINGS[system]['acd_ratio_range'][0], 
                            APHERESIS_SETTINGS[system]['acd_ratio_range'][1], 
                            acd_default)
    
    # Calculate results
    results = calculate_ecp_uv(tlc, lymph_percent, hct, system, lamp_power, target_dose,
                              use_hood, custom_distance, bag_type, flow_rate, plasma_removal, acd_ratio)
    
    # Display Results
    st.subheader(f"ECP Protocol for Hct {hct}% ({system})")
    
    # System Efficiency Panel
    eff_color = "red" if results['hct_efficiency'] < 0.85 else "green"
    st.markdown(f"""
    <div style="background-color:#f0f2f6;padding:10px;border-radius:5px;margin-bottom:20px">
        <h4 style="color:{eff_color}">System Efficiency: {results['hct_efficiency']:.2f} (1.0 = ideal at 40% Hct)</h4>
        <p>Hematocrit impact: <b>{APHERESIS_SETTINGS[system]['hct_impact']*100:.0f}%</b> sensitivity | 
        RBC contamination base: <b>{APHERESIS_SETTINGS[system]['rbc_base_contam']} ×10⁹</b></p>
    </div>
    """, unsafe_allow_html=True)
    
    # Results Columns
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("MNC Concentration", f"{results['mnc_conc']:.1f} ×10⁶/mL")
        st.metric("T-cell Purity", f"{results['purity_factor']:.2f}")
        
    with col2:
        st.metric("RBC Contamination", f"{results['rbc_contam']:.1f} ×10⁹", 
                delta=f"{(results['rbc_contam']-APHERESIS_SETTINGS[system]['rbc_base_contam']):.1f} vs baseline",
                delta_color="inverse")
        st.metric(f"{uv_type} Transmission", f"{results['transmission']*100:.1f}%")
        
    with col3:
        st.metric("Effective Dose", f"{results['effective_dose']:.2f} J/cm²")
        st.metric("Treatment Time", f"{results['exp_time']:.1f} min")
        st.metric("Distance", f"{results['distance']} cm")
    
    # Visualization
    st.subheader("Treatment Response Analysis")
    col1, col2 = st.columns(2)
    
    with col1:
        fig1, ax1 = plt.subplots(figsize=(8, 5))
        doses = np.linspace(0, target_dose*1.5, 100)
        ax1.plot(doses, 100*np.exp(-1.0*doses*results['transmission']*results['purity_factor']), 
               'r-', label='T-cells')
        ax1.plot(doses, 100*np.exp(-0.25*doses*results['transmission']), 
               'b-', label='Monocytes')
        ax1.axvline(results['effective_dose'], color='k', linestyle='--', label='Selected Dose')
        ax1.set_title(f'{uv_type} Dose-Response Curve')
        ax1.set_xlabel(f'{uv_type} Dose (J/cm²)')
        ax1.set_ylabel('Viability (%)')
        ax1.legend()
        ax1.grid(alpha=0.3)
        st.pyplot(fig1)
    
    with col2:
        fig2, ax2 = plt.subplots(figsize=(8, 5))
        times = np.linspace(0, max(results['exp_time']*2, 90), 100)
        time_doses = (results['effective_intensity']/1000) * (times * 60)
        ax2.plot(times, 100*np.exp(-1.0*time_doses*results['purity_factor']), 
               'r-', label='T-cells')
        ax2.plot(times, 100*np.exp(-0.25*time_doses), 
               'b-', label='Monocytes')
        ax2.axvline(results['exp_time'], color='k', linestyle='--', label='Estimated Time')
        ax2.set_title(f'{uv_type} Time-Response Curve')
        ax2.set_xlabel('Time (minutes)')
        ax2.set_ylabel('Viability (%)')
        ax2.legend()
        ax2.grid(alpha=0.3)
        st.pyplot(fig2)
    
    # Clinical Protocol
    st.subheader("Optimized Clinical Protocol")
    
    if hct > 45:
        st.warning("""
        **High Hematocrit Protocol Adjustments:**
        1. Reduce flow rate by 10-20% from standard
        2. Increase plasma removal by 5-10%
        3. Consider higher ACD ratio (1:14-1:16)
        """)
    
    st.markdown(f"""
    **For Hct {hct}% with {system}:**
    - **Flow Rate:** {flow_rate} mL/min ({'reduce' if hct > 45 else 'maintain'} for current Hct)
    - **ACD Ratio:** 1:{acd_ratio} ({'increase' if hct > 45 else 'standard'} anticoagulation)
    - **Expected RBC Contamination:** {results['rbc_contam']:.1f} ×10⁹ 
      ({(results['rbc_contam']/APHERESIS_SETTINGS[system]['rbc_base_contam']-1)*100:.0f}% {'increase' if hct > 40 else 'decrease'} from 40% Hct baseline)
    
    **UV Parameters:**
    - Maintain dose at {results['effective_dose']:.2f} ± 0.5 J/cm²
    - Treatment time: {results['exp_time']:.1f} minutes
    - Source distance: {results['distance']} cm
    """)

if __name__ == "__main__":
    main()
