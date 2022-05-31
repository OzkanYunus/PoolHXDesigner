#University of Bilkent
#ME430: HEAT EXCHANGERS § DESIGN
#SHELL AND TUBE HEAT EXCHANGER DESIGN FOR OLYMPIC POOLS
#Instructor: Dr. Barbaros ÇETİN
#Prepared by Yunus Ozkan and Bahadır Bilir
#Please read the instructions before use it.
#For any questions, please raise an issue on git or reach us by using e-mail, yunuszkn01@gmail.com

#All formulas taken from Heat Exchangers Selection, Rating and Thermal Design by Sadık Kakaç ,Hongtan Liu, Anchasa Pramuanjaroekenji

# Assume the tube length and outer diameter according to TEMA specifications
# 0.0666 < (Shell diameter/Tube length) < 0.2 (1)
# 1.25 < Pitch/Outer diameter of tube < 1.5 (2)


import pandas as pd
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import messagebox
import math



# This is given data for water

df_water = pd.read_csv("Water.csv", sep=",", header=None,
                       names=["Density", "Cp", "Viscosity", "K", "Pr", "Alpha"])


# This functions lineerly interpolate water at given temperature in between 0 to 300 Celcius Degree.
def waterLinterpolation(x):
    if x % 5 != 0 and x >= 0 and x <= 300:
        x_place = (x // 5) * 5
        x_remain = (x % 5)
        return ((df_water.loc[x_place + 5] * (x_remain)) + (df_water.loc[x_place] * (5 - x_remain))) / 5
    elif x % 5 == 0 and x >= 0 and x <= 300:
        return df_water.loc[x]
    else:
        tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a value in between 0-300 Celsius")


# These functions get values from buttons.

def buttonfunction1():
    global TColdIn
    try:
        a = float(Entry1.get())
        TColdIn = float(Entry1.get())
        if a <= 300 and a >= 0:
            return TColdIn
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a value in between 0-300 Celsius")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction2():
    global TColdOut
    try:
        a = float(Entry2.get())
        TColdOut = float(Entry2.get())
        if a <= 300 and a >= 0:
            return TColdOut
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a value in between 0-300 Celsius")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction3():
    global THotIn
    try:
        a = float(Entry3.get())
        THotIn = float(Entry3.get())
        if a <= 300 and a >= 0:
            return THotIn
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a value in between 0-300 Celsius")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def runprocess1():
    global TColdOut, TColdIn, Tbulk, Cold_properties
    try:
        Tbulk = (TColdOut + TColdIn) / 2

        Cold_properties = waterLinterpolation(Tbulk)
        messagebox.showinfo("Cold Water Bulk Properties", str(Cold_properties))


    except NameError:
        tk.messagebox.showwarning(title="UnknownValues", message="Please Provide All Required Data")
        Cold_properties = waterLinterpolation(Tbulk)
        messagebox.showinfo("Cold Water Bulk Properties", str(Cold_properties))


def buttonfunction5():
    global mdotcold
    try:
        a = float(Entry5.get())
        mdotcold = float(Entry5.get())
        if a <= 200 and a >= 0:
            return mdotcold
        else:
            tk.messagebox.showwarning(title="Invalid Interval",
                                      message="Please enter a value in between 0-250 kg per sec")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction6():
    global mdothot
    try:
        a = float(Entry6.get())
        mdothot = float(Entry6.get())
        if a <= 250 and a >= 0:
            return mdothot
        else:
            tk.messagebox.showwarning(title="Invalid Interval",
                                      message="Please enter a value in between 0-250 kg per sec")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


# 5 times iteration is used to find THotOut and its properties:

def runprocess2():
    global TColdOut, TColdIn, Tbulk, Cold_properties, THotIn, THotOut, mdothot, mdotcold, Qdot, Hot_properties, TbulkHot

    try:
        Qdot = mdotcold * Cold_properties.Cp * (TColdOut - TColdIn)
        Cp_hot = waterLinterpolation(THotIn).Cp  # Inital assumption is for ThotIn
        for i in range(5):
            THotOut = THotIn - (Qdot / (Cp_hot * mdothot))
            TbulkHot = (THotOut + THotIn) / 2
            Cp_hot = waterLinterpolation(TbulkHot).Cp
        Hot_properties = waterLinterpolation(TbulkHot)


        txtrunprocess2 = "Shell side hot water properties at {Tbulkhottemp:.2f} °C:" \
                         "\nHot water properties: {HotProperties}" \
                         "\nQ (Heat transfer): {HeatTransfer:.2f} kW " \
                         "\nHot fluid exit temperature: {Thotouttemp:.2f} °C".format(Tbulkhottemp = TbulkHot, HotProperties = Hot_properties, HeatTransfer = Qdot/1000, Thotouttemp = THotOut)
        messagebox.showinfo("Iterations Results", txtrunprocess2)

        #messagebox.showinfo("Hot Water Bulk Properties", "Iteration Results: " +
                            #"\n Hot Fluid Bulk Properties at:" + str(TbulkHot) + "Celsius\n" + str(Hot_properties) +
                            #"\n Q:" + str(Qdot) + "W" +"\n TShellout" + str(THotOut)) + "W"
    except NameError:
        tk.messagebox.showwarning(title="UnknownValues", message="Please Provide All Required Data")


def buttonfunction8():
    global d_out, r_out
    try:
        a = float(Entry8.get())
        d_out = float(Entry8.get()) / 1000
        r_out = d_out / 2

        if a <= 250 and a >= 0:
            return r_out
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a logical tube diameter value")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction9():
    global d_in, r_in
    try:
        a = float(Entry9.get())
        d_in = float(Entry9.get()) / 1000
        r_in = d_in / 2
        if a <= 250 and a >= 0:
            return r_in
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a logical tube diameter value")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction10():
    global Length
    try:
        a = float(Entry10.get())
        Length = float(Entry10.get())

        if a <= 35 and a >= 0:
            return Length
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a logical tube length in meter")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction11():
    global TColdInd
    try:
        a = float(Entry11.get())
        TColdInd = float(Entry11.get())
        if a <= 300 and a >= 0:
            return TColdInd
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a value in between 0-300 Celsius")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction12():
    global TColdOutd
    try:
        a = float(Entry12.get())
        TColdOutd = float(Entry12.get())
        if a <= 300 and a >= 0:
            return TColdOutd
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a value in between 0-300 Celsius")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction13():
    global THotInd
    try:
        a = float(Entry13.get())
        THotInd = float(Entry13.get())
        if a <= 300 and a >= 0:
            return THotInd
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a value in between 0-300 Celsius")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction14():
    global THotOutd
    try:
        a = float(Entry14.get())
        THotOutd = float(Entry14.get())
        if a <= 300 and a >= 0:
            return THotOutd
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a value in between 0-300 Celsius")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction15():
    global mdotcoldd
    try:
        a = float(Entry15.get())
        mdotcoldd = float(Entry15.get())
        if a <= 200 and a >= 0:
            return mdotcoldd
        else:
            tk.messagebox.showwarning(title="Invalid Interval",
                                      message="Please enter a value in between 0-250 kg per sec")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction16():
    global mdothotd
    try:
        a = float(Entry16.get())
        mdothotd = float(Entry16.get())
        if a <= 250 and a >= 0:
            return mdothotd
        else:
            tk.messagebox.showwarning(title="Invalid Interval",
                                      message="Please enter a value in between 0-250 kg per sec")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction17():
    global did
    try:
        a = float(Entry17.get())
        did = float(Entry17.get()) / 1000

        if a <= 250 and a >= 0:
            return did
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a logical tube diameter value")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction18():
    global dod
    try:
        a = float(Entry18.get())
        dod = float(Entry18.get()) / 1000

        if a <= 250 and a >= 0:
            return dod
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a logical tube diameter value")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction19():
    global Lengthd
    try:
        a = float(Entry19.get())
        Lengthd = float(Entry19.get())

        if a <= 35 and a >= 0:
            return Lengthd
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a logical tube length in meter")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction20():
    global bafflespacingd
    try:
        a = float(Entry20.get())
        bafflespacingd = float(Entry20.get())

        if a <= 35 and a >= 0:
            return bafflespacingd
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a logical length in meter")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction21():
    global Dsd
    try:
        a = float(Entry21.get())
        Dsd = float(Entry21.get())
        if a <= 300 and a >= 0:
            return Dsd
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a relevant value")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction22():
    global Ntd
    try:
        a = float(Entry22.get())
        Ntd = float(Entry22.get())
        if a <= 10000 and a >= 0:
            return Ntd
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a relevant value")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def buttonfunction23():
    global Tubesidepressurefd
    try:
        a = float(Entry23.get())
        Tubesidepressurefd = float(Entry23.get())
        if a <= 10000 and a >= 0:
            return Tubesidepressurefd
        else:
            tk.messagebox.showwarning(title="Invalid Interval", message="Please enter a relevant value")
    except ValueError:
        tk.messagebox.showwarning(title="ValueError", message="Please Provide Integer")


def runprocess3():
    global R_ft, h_hot_pre, h_cold_pre, r_in, r_out, Uf, Uc, k_material, lmtdf, lmtd, TColdOut, TColdIn, THotIn, THotOut, Qdot, Af, Ac, Nt, CTP, CL, Ds, PR, Length

    try:
        if choose_material.get() == "Titanium Alloy":
            R_ft = 0.000050
            k_material = 25
        elif choose_material.get() == "Stainless Steel":
            R_ft = 0.000090
            k_material = 15
        else:
            tk.messagebox.showwarning(title="UnknownValues", message="Please Provide All Required Data")

        if choose_pass.get() == "One Pass":
            CTP = 0.93

        elif choose_pass.get() == "Two Pass":
            CTP = 0.90

        else:
            tk.messagebox.showwarning(title="UnknownValues", message="Please Provide All Required Data")

        if var.get() == 1:
            CL = 1.0
        elif var.get() == 2:
            CL = 0.87
        else:
            tk.messagebox.showwarning(title="UnknownValues", message="Please Provide All Required Data")

        PR = float(pitch_ratio.get())

        h_hot_pre = float(htube_it.get())  # in W/(m^2.K)
        h_cold_pre = float(hshell_it.get())  # in W/(m^2.K)

        Uf = ((1 / h_hot_pre) + (r_out / r_in) * (1 / h_cold_pre) + R_ft + r_out * math.log(
            r_out / r_in) / k_material) ** (-1)
        Uc = ((1 / h_hot_pre) + (r_out / r_in) * (1 / h_cold_pre) + r_out * math.log(r_out / r_in) / k_material) ** (-1)
        try:
            lmtdf = ((THotIn - TColdOut) - (THotOut - TColdIn)) / (math.log((THotIn - TColdOut) / (THotOut - TColdIn)))
            lmtd = lmtdf * 0.9
            Af = Qdot / (Uf * lmtd)
            Ac = Qdot / (Uc * lmtd)
            Ds = 0.637 * math.sqrt(CL / CTP) * math.sqrt(Af * PR * PR * r_out * 2 / Length)
            Nt = (0.785 * ((CTP / CL) * Ds * Ds) / (PR * PR * r_out * r_out * 4))
        except RuntimeWarning:
            tk.messagebox.showwarning(title="UnknownValues",
                                      message="Since Delta T is very close to zero or zero, Error is occurred during" +
                                              "Calculations")
        txtrunprocess3 = "Suggested shell diameter: {SuggestShell:.3f} m " \
                         "\nRequired tube length: {TubeLength:.2f} m " \
                         "\nSuggested number of tubes: {NumberofTubes:.0f} " \
                         "\nGiven tube diameter OD:{ODDTube:.3f} m, ID: {IDDTube:.3f} " \
                         "\nSuggested baffle spacing for 25% cut: {suggestedbuffle:.3f} m" \
                         "\nSurface OverDesign: {overdesign:.3f}" \
                         "\nPitch Ratio: {PitchR:.2f}" \
                         "\nHot water shell exit temperature {HotexitTemperature:.2f} °C".format(SuggestShell = Ds, TubeLength = Length, NumberofTubes = Nt, ODDTube = r_out*2 ,IDDTube = r_in*2, suggestedbuffle= Ds*0.6, overdesign = Af/Ac, PitchR = PR, HotexitTemperature = THotOut)

        messagebox.showinfo("Preliminary Design Results", txtrunprocess3)

        # messagebox.showinfo("Preliminary Design Results", "Suggested Shell Diameter in m:   " + str(Ds) +
        #                     "\nRequired Tube Length in m:   " + str(Length) +
        #                     "\nMax Number of Tubes:   " + str(Nt) +
        #                     "\nGiven Tube Diameter OD in m: " + str(r_out * 2) + "ID:" + str(r_in * 2) +
        #                     "\nSuggested Baffle Spacing for 25% cut in m: " + str(0.6 * Ds) +
        #                     "\nSurface Over Design:   " + str(Af / Ac) +
        #                     "\n Pitch Ratio:   " + str(PR) +
        #                     "\nShell out Temperature in Celsius:   " + str(THotOut))

    except NameError:
        tk.messagebox.showwarning(title="UnknownValues", message="Please Provide All Required Data")


def runprocess4():
    global Atpd, Ded, Asd, Cd, Gsd, Re_sd, h_td, h_sd, TColdOutd, TColdInd, TbulkColdd, TbulkHotd, pitch_sized, mdothotd, mdotcoldd

    # Calculating tube side cold water heat transfer coefficient:

    if choose_passd.get() == "One Pass":
        NoPd = 1
    else:  # if choose_passd.get() == "Two Pass": if you want to add more pass option
        NoPd = 2

    pitch_sized = float(pitch_ratiod.get()) * dod

    TbulkColdd = (TColdInd + TColdOutd) / 2

    Cold_propertiesd = waterLinterpolation(TbulkColdd)

    Atpd = (math.pi * (did) ** 2) / 4 * Ntd / NoPd

    Umd = mdotcoldd / (Cold_propertiesd.Density * Atpd)

    Re_td = (Cold_propertiesd.Density * Umd * did) / Cold_propertiesd.Viscosity

    # Gnielinski's correlation for inside the tube, cold flow.

    if Re_td > 10 ** 4:
        f_td = (1.58 * math.log(Re_td) - 3.28) ** -2

        Nu_td = ((f_td / 2) * (Re_td - 1000) * Cold_propertiesd.Pr) / (
        (1 + 12.7 * math.sqrt(f_td / 2) * ((Cold_propertiesd.Pr) ** (2 / 3) - 1)))

        h_td = (Nu_td * Cold_propertiesd.K) / did



    else:

        f_td = (1.58 * math.log(Re_td) - 3.28) ** -2

        Nu_td = ((f_td / 2) * (Re_td - 1000) * Cold_propertiesd.Pr) / (
        (1 + 12.7 * math.sqrt(f_td / 2) * ((Cold_propertiesd.Pr) ** (2 / 3) - 1)))

        h_td = (Nu_td * Cold_propertiesd.K) / did

        tk.messagebox.showwarning(title="Low Reynolds Number",
                                  message="This program is using Gnielinski's correlation, Re number should be greater than 10^4 for cold flow." +
                                          "\n Displayed results will be wrong."
                                          "\n Try to increase the flow rate or try smaller tubes.")

    # Shell side hot water heat transfer coefficient

    TbulkHotd = (THotOutd + THotInd) / 2

    Hot_propertiesd = waterLinterpolation(TbulkHotd)

    Cd = pitch_sized - dod

    if choose_pitchLayout.get() == "Square Pitch":
        Ded = 4 * ((pitch_sized ** 2) - (math.pi * (dod ** 2) / 4)) / (math.pi * dod)

    elif choose_pitchLayout.get() == "Triangular Pitch":
        Ded = 4 * (((pitch_sized ** 2) * math.sqrt(3) / 4) - (math.pi * (dod ** 2) / 8)) / (math.pi * dod / 2)


    Asd = Dsd * Cd * bafflespacingd / pitch_sized
    Gsd = mdothotd / Asd

    Re_sd = Gsd * Ded / Hot_propertiesd.Viscosity

    if Re_sd > 10 ** 6 or Re_sd < 2 * 10 ** 3:
        tk.messagebox.showwarning(title="Low Reynolds Number",
                                  message="This program is using Kernn Method, Shell side Re number should be in between 2x10^3 - 10^6" +
                                          "\n Displayed results will be wrong.")

    T_wd = 1 / 2 * ((TColdInd + TColdOutd) / 2 + (THotInd + THotOutd) / 2)
    wall_property_sd = waterLinterpolation(T_wd)
    h_sd = Hot_propertiesd.K * 0.36 * (((Ded * Gsd / Hot_propertiesd.Viscosity) ** 0.55) * (
                (Hot_propertiesd.Cp * Hot_propertiesd.Viscosity / Hot_propertiesd.K) ** (1 / 3)) * ((
                                                                                                                Hot_propertiesd.Viscosity / wall_property_sd.Viscosity) ** 0.14)) / Ded

    # Use these results to iterate h values in runprocess 3 until difference between h values is insignificant.

    txtrunprocess4 = "Tube side heat transfer coefficient: {htubeside: .2f} W/m^2*K" \
                     "\nShell side heat transfer coefficient: {hshellside: .2f} W/m^2*K" \
                     "\nTube side reynolds number: {Rtubeside: .2f} " \
                     "\nShell side reynolds number: {Rshellside: .2f}".format(htubeside=h_td, hshellside=h_sd, Rtubeside=Re_td, Rshellside=Re_sd)

    messagebox.showinfo("Results", txtrunprocess4)

    # messagebox.showinfo("Results", "h tube side:   " + str(h_td) + "W/m^2*K" + "h shell side:   "+ str(h_sd) + "W/m^2*K"
    #                     "\n, Tube side Reynolds number:   " + str(Re_td) + "Shell side Reynolds number" + str(Re_sd))


def runprocess5():
    global Atpd, Ded, Asd, Cd, Gsd, Re_sd, h_td, h_sd, TColdOutd, TColdInd, TbulkColdd, TbulkHotd, pitch_sized, mdothotd, mdotcoldd

    # Calculating tube side cold water heat transfer coefficient:

    if choose_passd.get() == "One Pass":
        NoPd = 1
    else:  # if choose_passd.get() == "Two Pass": if you want to add more pass option
        NoPd = 2

    pitch_sized = float(pitch_ratiod.get()) * dod

    TbulkColdd = (TColdInd + TColdOutd) / 2

    Cold_propertiesd = waterLinterpolation(TbulkColdd)

    Atpd = (math.pi * (did) ** 2) / 4 * Ntd / NoPd

    Umd = mdotcoldd / (Cold_propertiesd.Density * Atpd)

    Re_td = (Cold_propertiesd.Density * Umd * did) / Cold_propertiesd.Viscosity

    # Gnielinski's correlation for inside the tube, cold flow.

    if Re_td > 10 ** 4:
        f_td = (1.58 * math.log(Re_td) - 3.28) ** -2

        Nu_td = ((f_td / 2) * (Re_td - 1000) * Cold_propertiesd.Pr) / (
        (1 + 12.7 * math.sqrt(f_td / 2) * ((Cold_propertiesd.Pr) ** (2 / 3) - 1)))

        h_td = (Nu_td * Cold_propertiesd.K) / did



    else:

        f_td = (1.58 * math.log(Re_td) - 3.28) ** -2

        Nu_td = ((f_td / 2) * (Re_td - 1000) * Cold_propertiesd.Pr) / (
        (1 + 12.7 * math.sqrt(f_td / 2) * ((Cold_propertiesd.Pr) ** (2 / 3) - 1)))

        h_td = (Nu_td * Cold_propertiesd.K) / did

        tk.messagebox.showwarning(title="Low Reynolds Number",
                                  message="This program is using Gnielinski's correlation, Re number should be greater than 10^4 for cold flow." +
                                          "\n Displayed results will be wrong."
                                          "\n Try to increase the flow rate or try smaller tubes.")

    # Shell side hot water heat transfer coefficient

    TbulkHotd = (THotOutd + THotInd) / 2

    Hot_propertiesd = waterLinterpolation(TbulkHotd)

    Cd = pitch_sized - dod

    if choose_pitchLayout.get() == "Square Pitch":
        Ded = 4 * ((pitch_sized ** 2) - (math.pi * (dod ** 2) / 4)) / (math.pi * dod)

    elif choose_pitchLayout.get() == "Triangular Pitch":
        Ded = 4 * (((pitch_sized ** 2) * math.sqrt(3) / 4) - (math.pi * (dod ** 2) / 8)) / (math.pi * dod / 2)


    Asd = Dsd * Cd * bafflespacingd / pitch_sized
    Gsd = mdothotd / Asd

    Re_sd = Gsd * Ded / Hot_propertiesd.Viscosity

    if Re_sd > 10 ** 6 or Re_sd < 2 * 10 ** 3:
        tk.messagebox.showwarning(title="Low Reynolds Number",
                                  message="This program is using Kernn Method, Shell side Re number should be in between 2x10^3 - 10^6" +
                                          "\n Displayed results will be wrong.")

    T_wd = 1 / 2 * ((TColdInd + TColdOutd) / 2 + (THotInd + THotOutd) / 2)
    wall_property_sd = waterLinterpolation(T_wd)
    h_sd = Hot_propertiesd.K * 0.36 * (((Ded * Gsd / Hot_propertiesd.Viscosity) ** 0.55) * (
                (Hot_propertiesd.Cp * Hot_propertiesd.Viscosity / Hot_propertiesd.K) ** (1 / 3)) * ((
                                                                                                                Hot_propertiesd.Viscosity / wall_property_sd.Viscosity) ** 0.14)) / Ded

    # If Iteration Completed checkbox is not checked:

    # If Iteration Completed checkbox is checked:
    # U calculations:
    if choose_materiald.get() == "Titanium Alloy":
        R_ftd = 0.000050
        k_materiald = 25
    elif choose_materiald.get() == "Stainless Steel":
        R_ftd = 0.000090
        k_materiald = 15
    else:
        tk.messagebox.showwarning(title="UnknownValues", message="Please Provide All Required Data")

    Ufd = ((1 / h_sd) + ((dod / 2) / (did / 2)) * (1 / h_td) + R_ftd + (dod / 2) * math.log(
        (dod / 2) / (did / 2)) / k_materiald) ** (-1)

    Ucd = ((1 / h_sd) + ((dod / 2) / (did / 2)) * (1 / h_td) + (dod / 2) * math.log(
        (dod / 2) / (did / 2)) / k_materiald) ** (-1)


    # Calculate pressure drop
    # Shellside pressure drop

    f_pressure_sd = math.exp(0.576 - 0.19 * math.log(Re_sd))

    flux_factor = (Hot_propertiesd.Viscosity / wall_property_sd.Viscosity) ** 0.14
    NumberofBaffles = (Lengthd / bafflespacingd) - 1

    shellsidePressurechange = (f_pressure_sd * (Gsd ** 2) * (NumberofBaffles + 1) * Dsd) / (
                Hot_propertiesd.Density * Ded * flux_factor)

    # Tubeside pressure drop

    tubesidePressurechange = (4 * Tubesidepressurefd * (NoPd) * (Gsd ** 2)) / (did * (Cold_propertiesd.Density))

    txtrunprocess5 = "Tube side heat transfer coefficient: {h_tubeside:.2f} (W/m^2*K)   Shell side heat transfer coefficient: {h_shellside:.2f} (W/m^2*K) " \
                     "\nTube side Reynolds Number: {R_tubeside:.2f}    Shell side Reynolds Number: {R_shellside:.2f} " \
                     "\nOverall heat transfer coefficient: {U_all:.2f} (W/m^2*K)   Over Surface Design: {OS:.2f} " \
                     "\n ube side pressure drop: {P_tubeside:.2f} kPa   Shell side pressure drop: {P_shellside:.2f} kPa".format(h_tubeside = h_td, h_shellside = h_sd, R_tubeside = Re_td, R_shellside = Re_sd, U_all = Ufd, OS = Ucd/Ufd, P_tubeside = tubesidePressurechange/1000, P_shellside = shellsidePressurechange/1000)

    ttk.Label(frame_right2, text = txtrunprocess5).pack(padx=3, pady=0, side="left")




    # ttk.Label(frame_right2, text="h tube side: " + str(h_td) + "(W/m^2*K)"+ "h shell side: " + str(h_sd) + "(W/m^2*K)"+
    #                     "\n, Tube side Reynolds number: " + str(Re_td) + "Shell side Reynolds number: " + str(Re_sd)+
    #                     "\n Overall heat transfer coefficient: "+  str(Ufd) + "(W/m^2*K)"  + "Over surface design:   "+ str(Ucd/Ufd)+
    #                     "\n Tube side pressure change    " + str(tubesidePressurechange/1000) + "kPa" + "Shell side pressure change   "+str(shellsidePressurechange/1000) +"kPa").pack(padx=3, pady=0, side="left")



# User Interface

# Top and mid frames are for preliminary design

root = tk.Tk()
root.title("Water to Water POOL HX DESIGNER")
canvas = tk.Canvas(root, height=1000, width=1650)
canvas.pack()


style = ttk.Style(root)
style.theme_use('clam')
style.configure('TButton',font=('Arial', 10, "bold"), foreground= "black", background= 'gainsboro')
style.configure('TLabel',font=('Arial', 10), foreground= "black", background= 'white' )
#style.configure('TEntry', foreground= "red", background= 'black' )
style.configure('TRadiobutton', foreground= "black", background= 'white' )




frame_top = tk.Frame(root, bg="#add8e6")
frame_top.place(relx=0.05, rely=0.05, relwidth=0.90, relheight=0.1)

frame_mid = tk.Frame(root, bg="#add8e6")
frame_mid.place(relx=0.05, rely=0.16, relwidth=0.90, relheight=0.1)

frame_mid2 = tk.Frame(root, bg="aqua")
frame_mid2.place(relx=0.05, rely=0.27, relwidth=0.90, relheight=0.05)

# Left frames for main process

frame_left = tk.Frame(root, bg="rosybrown")
frame_left.place(relx=0.05, rely=0.35, relwidth=0.30, relheight=0.10)

frame_left2 = tk.Frame(root, bg="rosybrown")
frame_left2.place(relx=0.05, rely=0.45, relwidth=0.30, relheight=0.05)

frame_left3 = tk.Frame(root, bg="rosybrown")
frame_left3.place(relx=0.05, rely=0.50, relwidth=0.30, relheight=0.05)

frame_left4 = tk.Frame(root, bg="rosybrown")
frame_left4.place(relx=0.05, rely=0.55, relwidth=0.30, relheight=0.05)

frame_left5 = tk.Frame(root, bg="rosybrown")
frame_left5.place(relx=0.05, rely=0.60, relwidth=0.30, relheight=0.05)

frame_left6 = tk.Frame(root, bg="rosybrown")
frame_left6.place(relx=0.05, rely=0.65, relwidth=0.30, relheight=0.05)

frame_left7 = tk.Frame(root, bg="rosybrown")
frame_left7.place(relx=0.05, rely=0.70, relwidth=0.30, relheight=0.05)

frame_left8 = tk.Frame(root, bg="rosybrown")
frame_left8.place(relx=0.05, rely=0.75, relwidth=0.30, relheight=0.05)

frame_left9 = tk.Frame(root, bg="rosybrown")
frame_left9.place(relx=0.05, rely=0.80, relwidth=0.30, relheight=0.05)

Main_Calculation = tk.Label(frame_left, bg= "rosybrown", text="Main Design Calculations", font="Verdana 10 bold").pack(padx=10, pady=25,
                                                                                                side="top")

frame_right = tk.Frame(root, bg="lightcoral")
frame_right.place(relx=0.40, rely=0.35, relwidth=0.55, relheight=0.10)

frame_right2 = tk.Frame(root, bg="lightcoral")
frame_right2.place(relx=0.40, rely=0.45, relwidth=0.55, relheight=0.40)

constraints = tk.Label(frame_top, bg="#add8e6", text="Fluid Property Calculations", font="Verdana 10 bold").pack(padx=10, pady=5,
                                                                                                side="top")
constraints2 = tk.Label(frame_mid, bg="#add8e6", text="Preliminary Design", font="Verdana 10 bold").pack(padx=10,
                                                                                                         pady=5,
                                                                                                         side="top")
ttk.Label(frame_top, text="Tube Inlet °C").pack(padx=5, pady=5, side="left")
Entry1 = tk.Entry(frame_top, width=7)
Entry1.pack(padx=1, pady=0, side="left")
ttk.Button(frame_top, text="Send", command=buttonfunction1).pack(padx=1, pady=0, side="left")
ttk.Label(frame_top, text="Tube Outlet °C").pack(padx=3, pady=0, side="left")

Entry2 = ttk.Entry(frame_top, width=7)
Entry2.pack(padx=1, pady=0, side="left")

ttk.Button(frame_top, text="Send", command=buttonfunction2).pack(padx=1, pady=0, side="left")
ttk.Button(frame_top, text="Tubeside properties", command=runprocess1).pack(padx=5, pady=0, side="left")
ttk.Label(frame_top, text="Shell Inlet °C").pack(padx=3, pady=0, side="left")
Entry3 = ttk.Entry(frame_top, width=7)
Entry3.pack(padx=1, pady=0, side="left")
ttk.Button(frame_top, text="Send", command=buttonfunction3).pack(padx=3, pady=0, side="left")
ttk.Label(frame_top, text="Tube (kg/s)").pack(padx=3, pady=0, side="left")
Entry5 = ttk.Entry(frame_top, width=7)
Entry5.pack(padx=1, pady=0, side="left")
ttk.Button(frame_top, text="Send", command=buttonfunction5).pack(padx=3, pady=0, side="left")
tk.Label(frame_top, text="Shell (kg/s)").pack(padx=3, pady=0, side="left")
Entry6 = ttk.Entry(frame_top, width=7)
Entry6.pack(padx=1, pady=0, side="left")
ttk.Button(frame_top, text="Send", command=buttonfunction6).pack(padx=3, pady=0, side="left")
ttk.Button(frame_top, text="Property Results", command=runprocess2).pack(padx=25, pady=0, side="left")

# Preliminary Design UI
ttk.Label(frame_mid,  text="Tube OD in mm").pack(padx=3, pady=0, side="left")
Entry8 = tk.Entry(frame_mid, width=7)
Entry8.pack(padx=1, pady=0, side="left")
ttk.Button(frame_mid, text="Send", command=buttonfunction8).pack(padx=3, pady=0, side="left")
ttk.Label(frame_mid,  text="Tube ID in mm").pack(padx=3, pady=0, side="left")
Entry9 = ttk.Entry(frame_mid, width=7)
Entry9.pack(padx=1, pady=0, side="left")
ttk.Button(frame_mid, text="Send", command=buttonfunction9).pack(padx=3, pady=0, side="left")
ttk.Label(frame_mid,  text="Tube Length in m").pack(padx=3, pady=0, side="left")
Entry10 = tk.Entry(frame_mid, width=7)
Entry10.pack(padx=1, pady=0, side="left")
ttk.Button(frame_mid, text="Send", command=buttonfunction10).pack(padx=3, pady=0, side="left")

choose_material = tk.StringVar(frame_mid)
choose_material.set("Set Material")
choose_material_menu = tk.OptionMenu(frame_mid, choose_material, "Titanium Alloy", "Stainless Steel")
choose_material_menu.pack(padx=3, pady=0, side="left")

choose_pass = tk.StringVar(frame_mid)
choose_pass.set("Set Number of Pass")
choose_pass_menu = tk.OptionMenu(frame_mid, choose_pass, "One Pass", "Two Pass")
choose_pass_menu.pack(padx=15, pady=0, side="top", anchor="nw")

var = tk.IntVar()
ttk.Radiobutton(frame_mid, text="90° or 45°", variable=var, value=1).pack(padx=3, pady=0, side="left")
ttk.Radiobutton(frame_mid, text="30° or 60°", variable=var, value=2).pack(padx=3, pady=0, side="left")

ttk.Label(frame_mid, text="Pitch Ratio:").pack(padx=3, pady=0, side="left")

pitch_ratio = tk.StringVar(frame_mid)
pitch_ratio.set("1.25")
choose_ratio_menu = tk.OptionMenu(frame_mid, pitch_ratio, "1.25", "1.30", "1.35", "1.40", "1.45", "1.50")
choose_ratio_menu.pack(padx=3, pady=0, side="left", anchor="nw")

htube_it = tk.StringVar(frame_mid2)
htube_it.set("5000")
ttk.Label(frame_mid2, text="h tube iteration:").pack(padx=3, pady=0, side="left")
htube_it_menu = tk.OptionMenu(frame_mid2, htube_it, "2000", "2100", "2200", "2300", "2400", "2500", "2600", "2700",
                               "2800", "2900","3000", "3100", "3200", "3300", "3400", "3500", "3600", "3700",
                              "3800", "3900",
                              "4000", "4100", "4200", "4300", "4400", "4500", "4600", "4700", "4800", "4900",
                              "5000", "5100", "5200", "5300", "5400", "5500", "5600", "5700", "5800", "5900",
                              "6000", "6100", "6200", "6300", "6400", "6500", "6600", "6700", "6800", "6900")
htube_it_menu.pack(padx=15, pady=0, side="left")

hshell_it = tk.StringVar(frame_mid2)
hshell_it.set("5000")
ttk.Label(frame_mid2, text="h shell iteration:").pack(padx=3, pady=0, side="left")
hshell_it_menu = tk.OptionMenu(frame_mid2, hshell_it, "2000", "2100", "2200", "2300", "2400", "2500", "2600", "2700",
                               "2800", "2900",
                               "3000", "3100", "3200", "3300", "3400", "3500", "3600", "3700", "3800", "3900",
                               "4000", "4100", "4200", "4300", "4400", "4500", "4600", "4700", "4800", "4900",
                               "5000", "5100", "5200", "5300", "5400", "5500", "5600", "5700", "5800", "5900",
                               "6000", "6100", "6200", "6300", "6400", "6500", "6600", "6700", "6800", "6900")

hshell_it_menu.pack(padx=15, pady=0, side="left")

ttk.Button(frame_mid2, text="Preliminary Design Results", command=runprocess3).pack(padx=10, pady=0, side="left")

ttk.Label(frame_mid2,
         text="Initial values are 5000. Change h, with respect to main part results and change main values with respect to"
              "\n preliminary design results. Press end of the iteration, when you find out difference in h values are less than 100").pack(
    padx=20, pady=0, side="left")

# Main UI
ttk.Label(frame_left, text="Tube Inlet °C").pack(padx=3, pady=0, side="left", anchor="nw")
Entry11 = ttk.Entry(frame_left, width=7)
Entry11.pack(padx=1, pady=0, side="left", anchor="nw")
ttk.Button(frame_left, text="Send", command=buttonfunction11).pack(padx=3, pady=0, side="left", anchor="nw")

ttk.Label(frame_left, text="Tube Outlet °C").pack(padx=3, pady=0, side="left", anchor="nw")
Entry12 = ttk.Entry(frame_left, width=7)
Entry12.pack(padx=1, pady=0, side="left", anchor="nw")
ttk.Button(frame_left, text="Send", command=buttonfunction12).pack(padx=3, pady=0, side="left", anchor="nw")

ttk.Label(frame_left2,  text="Shell Inlet °C").pack(padx=3, pady=0, side="left")
Entry13 = ttk.Entry(frame_left2, width=7)
Entry13.pack(padx=1, pady=0, side="left")
ttk.Button(frame_left2, text="Send", command=buttonfunction13).pack(padx=3, pady=0, side="left")

ttk.Label(frame_left2, text="Shell Outlet °C").pack(padx=3, pady=0, side="left")
Entry14 = ttk.Entry(frame_left2, width=7)
Entry14.pack(padx=1, pady=0, side="left")
ttk.Button(frame_left2, text="Send", command=buttonfunction14).pack(padx=3, pady=0, side="left")

ttk.Label(frame_left3, text="Tube (kg/s)").pack(padx=3, pady=0, side="left", anchor="nw")
Entry15 = ttk.Entry(frame_left3, width=7)
Entry15.pack(padx=3, pady=0, side="left", anchor="nw")
ttk.Button(frame_left3, text="Send", command=buttonfunction15).pack(padx=3, pady=0, side="left", anchor="nw")

ttk.Label(frame_left3, text="Shell (kg/s)").pack(padx=10, pady=0, side="left", anchor="nw")
Entry16 = ttk.Entry(frame_left3, width=7)
Entry16.pack(padx=1, pady=0, side="left", anchor="nw")
ttk.Button(frame_left3, text="Send", command=buttonfunction16).pack(padx=3, pady=0, side="left", anchor="nw")

# Entry 18-17 re-placed for better screening

ttk.Label(frame_left4, text="Tube OD in mm").pack(padx=0, pady=0, side="left", anchor="nw")
Entry18 = ttk.Entry(frame_left4, width=7)
Entry18.pack(padx=1, pady=0, side="left", anchor="nw")
ttk.Button(frame_left4, text="Send", command=buttonfunction18).pack(padx=3, pady=0, side="left", anchor="nw")

ttk.Label(frame_left4, text="Tube ID in mm").pack(padx=0, pady=0, side="left", anchor="nw")
Entry17 = ttk.Entry(frame_left4, width=7)
Entry17.pack(padx=3, pady=0, side="left", anchor="nw")
ttk.Button(frame_left4, text="Send", command=buttonfunction17).pack(padx=3, pady=0, side="left", anchor="nw")

ttk.Label(frame_left5, text="Tube length in m").pack(padx=0, pady=0, side="left", anchor="nw")
Entry19 = ttk.Entry(frame_left5, width=7)
Entry19.pack(padx=3, pady=0, side="left", anchor="nw")
ttk.Button(frame_left5, text="Send", command=buttonfunction19).pack(padx=3, pady=0, side="left", anchor="nw")

ttk.Label(frame_left5, text="Baffle Spacing in m").pack(padx=0, pady=0, side="left", anchor="nw")
Entry20 = ttk.Entry(frame_left5, width=7)
Entry20.pack(padx=1, pady=0, side="left", anchor="nw")
ttk.Button(frame_left5, text="Send", command=buttonfunction20).pack(padx=3, pady=0, side="left", anchor="nw")

ttk.Label(frame_left6, text="Shell Diameter").pack(padx=0, pady=0, side="left", anchor="nw")
Entry21 = ttk.Entry(frame_left6, width=7)
Entry21.pack(padx=3, pady=0, side="left", anchor="nw")
ttk.Button(frame_left6, text="Send", command=buttonfunction21).pack(padx=3, pady=0, side="left", anchor="nw")

ttk.Label(frame_left6, text="Number of Tubes").pack(padx=0, pady=0, side="left", anchor="nw")
Entry22 = ttk.Entry(frame_left6, width=7)
Entry22.pack(padx=1, pady=0, side="left", anchor="nw")
ttk.Button(frame_left6, text="Send", command=buttonfunction22).pack(padx=3, pady=0, side="left", anchor="nw")

choose_pitchLayout = tk.StringVar(frame_left7)
choose_pitchLayout.set("Square Pitch")
ttk.Label(frame_left7, text="Pitch Type:").pack(padx=5, pady=0, side='left', anchor="nw")
choose_pitch_menu = tk.OptionMenu(frame_left7, choose_pitchLayout, "Square Pitch", "Triangular Pitch")
choose_pitch_menu.pack(pady=0, side='left', anchor="nw")

choose_materiald = tk.StringVar(frame_left7)
choose_materiald.set("Titanium Alloy")
ttk.Label(frame_left7, text="Material:").pack(padx=5, pady=0, side='left', anchor="nw")
choose_materiald_menu = tk.OptionMenu(frame_left7, choose_materiald, "Titanium Alloy", "Stainless Steel")
choose_materiald_menu.pack(padx=0, pady=0, side="left", anchor="nw")

pitch_ratiod = tk.StringVar(frame_left8)
pitch_ratiod.set("1.25")
ttk.Label(frame_left8, text="Pitch Ratio").pack(padx=5, pady=0, side='left', anchor="nw")
choose_ratiod_menu = tk.OptionMenu(frame_left8, pitch_ratiod, "1.25", "1.30", "1.35", "1.40", "1.45", "1.50")
choose_ratiod_menu.pack(padx=3, pady=0, side="left", anchor="nw")

choose_passd = tk.StringVar(frame_left8)
choose_passd.set("One Pass")
ttk.Label(frame_left8, text="Number of Pass").pack(padx=5, pady=0, side='left', anchor="nw")
choose_passd_menu = tk.OptionMenu(frame_left8, choose_passd, "One Pass", "Two Pass")
choose_passd_menu.pack(padx=3, pady=0, side="left", anchor="nw")

ttk.Button(frame_left9, text="Display h values ", command=runprocess4).pack(padx=150, pady=0, side="top", anchor="nw")

ttk.Label(frame_right, text="Tube Side Pressure Drop Correction Factor:").pack(padx=0, pady=0, side="left", anchor="nw")
Entry23 = ttk.Entry(frame_right, width=7)
Entry23.pack(padx=1, pady=0, side="left", anchor="nw")
ttk.Button(frame_right, text="Send", command=buttonfunction23).pack(padx=3, pady=0, side="left", anchor="nw")

ttk.Button(frame_right, text="Display All Results ", command=runprocess5).pack(padx=150, pady=0, side="left", anchor="nw")

root.mainloop()
