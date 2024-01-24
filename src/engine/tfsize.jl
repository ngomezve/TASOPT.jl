"""
    tfsize!(gee, M0, T0, p0, a0, M2, M25,
      Feng, Phiinl, Kinl, iBLIc,
      BPR, pif, pilc, pihc,
      pid, pib, pifn, pitn,
      Ttf, ifuel, etab,
      epf0, eplc0, ephc0, epht0, eplt0,
      pifK, epfK,
      mofft, Pofft,
      Tt9, pt9, Tt4,
      epsl, epsh,
      icool,
      Mtexit, dTstrk, StA, efilm, tfilm,
      M4a, ruc,
      ncrowx, ncrow,
      epsrow, Tmrow, 
      Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
      Δp_PreC, Δp_InterC, Δp_Regen)


Turbofan performance and sizing routine.
      
Calculation procedure follows that of Kerrebrock, but the usual gas property formulas are replaced by function calls, which can therefore implement more general gas models.  
In addition, a turbine cooling model is added.
      
The gas routines reside in the following source files:
    gascalc.f  Routines for various processes (compressor, turbine, combustor, etc)
    gasfun.f   Routines for computing cp[T], h[t], sigma[T], R, called by the routines in gascalc.f
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
    - `gee`:     gravity acceleration
    - `M0`:      freestream Mach
    - `T0`:      freestream temperature  [K]
    - `p0`:      freestream pressure  [Pa]
    - `M2`:      fan-face Mach number
    - `M25`:     HPC-face Mach number
    - `Feng`:    required net thrust  (PK_inl+PK_out-Phi_jet)/u0  =  sum( mdot u)
    - `Phiinl`:  inlet ingested dissipation
    - `iBLIc`:   0=core in clear flow, 1=core sees Phiinl
    - `BPR`:     bypass ratio  = mdot_fan/mdot_core
    - `pif`:     fan      pressure ratio  ( = pt7 /pt2)
    - `pilc`:    LP comp  pressure ratio  ( = pt25/pt2)
    - `pihc`:    HP comp  pressure ratio  ( = pt3 /pt25)
    - `pid`:     diffuser pressure ratio  ( = pt2 /pt0)
    - `pib`:     burner   pressure ratio  ( = pt4 /pt3)
    - `pifn`:    fan     nozzle pressure ratio  ( = pt7/pt2.1)
    - `pitn`:    turbine nozzle pressure ratio  ( = pt5/pt4.9)
    - `Ttf`:     fuel temperature entering combustor
    - `ifuel`:   fuel index, see function gasfun (in gasfun.f)
    - `etab`:    combustor efficiency (fraction of fuel burned)
    - `epf0`:    fan max polytropic efficiency
    - `eplc0`:   LPC max polytropic efficiency
    - `ephc0`:   HPC max polytropic efficiency
    - `epht0`:   HPT max polytropic efficiency
    - `eplt0`:   LPT max polytropic efficiency
    - `pifK`:    fan efficiency FPR offset:    epolf = epf0 + epfK*(pif-pifK)
    - `epfK`:    fan efficiency pif derivative
      
    - `mofft`:    mass flow offtake at LPC discharge station 2.5
    - `Pofft`:    low spool power offtake
    - `Tt9`:     offtake air discharge total temperature
    - `pt9`:     offtake air discharge total pressure
    - `epsl`:    low  spool power loss fraction
    - `epsh`:    high spool power loss fraction
      
    - `icool   turbine cooling flag
               0 = no cooling, ignore all cooling parameters below
               1 = usual cooling, using passed-in fcool
               2 = usual cooling, but set (and return) fcool from Tmetal
    - `Mtexit`:   turbine blade-row exit Mach, for setting temperature drops
    - `dTstrk`:   hot-streak temperature delta {K}, used only if icool=2
    - `StA`:      area-weighted Stanton number    , used only if icool=2
    - `M4a`:      effective Mach at cooling-flow outlet (start of mixing)
    - `ruc`:      cooling-flow outlet velocity ratio, u/ue
    - `ncrowx`:      dimension of epsrow array
    - `ncrow`:       number of blade rows requiring cooling
    - `epsrow(.)`:   input specified  cooling-flow bypass ratio if icool=1
                     output resulting cooling-flow bypass ratio if icool=2
    - `Tmrow(.)`:    input specified  metal temperature  [K]    if icool=2
                     output resulting metal temperature  [K]    if icool=1

      **Outputs:**
    - `epsrow(.)`:   see above
    - `Tmrow(.)`:    see above
    - `TSFC`:    thrust specific fuel consumption = mdot_fuel g / F   [1/s]
    - `Fsp`:     specific thrust  = F / (mdot u0) = F / ((1+BPR) mdot_core u0)
    - `hfuel`:   fuel heating value   [J / kg K]
    - `ff`:      fuel mass flow fraction  =  mdot_fuel / mdot_core
    - `mcore`:   core mass flow = mdot_core  [kg/s]
    - `A2`:      fan-face area [m^2]
    - `A25`:     HPC-face area [m^2]
    - `A5`:      core nozzle area [m^2]
    - `A7`:      fan  nozzle area [m^2]
    - `A6`:      core plume  area [m^2]
    - `A8`:      fan  plume  area [m^2]
    - `Tt?`:     total temperature
    - `ht?`:     total complete enthalpy (includes heat of formation)
    - `pt?`:     total pressure
    - `cpt?`:    specific heat at stagnation temperature  (= dh/dT)
    - `Rt?`:     gas constant  at stagnation conditions
    - `T?`:      static temperature
    - `u?`:      velocity
    - `epf`:     fan polytropic efficiency
    - `eplc`:    LPC polytropic efficiency
    - `ephc`:    HPC polytropic efficiency
    - `epht`:    HPT polytropic efficiency
    - `eplt`:    LPT polytropic efficiency
    - `etaf`:    fan overall efficiency
    - `etalc`:   LPC overall efficiency
    - `etahc`:   HPC overall efficiency
    - `etaht`:   HPT overall efficiency
    - `etalt`:   LPT overall efficiency
    - `Lconv`:   T if convergence was successful, F otherwise

    The "?" symbol denotes the station index:
      0  freestream
      18 fan face outside of casing BLs
      19 fan face over LPC portion
      2  fan face over fan portion
      21 fan exit
      25 LPC exit, HPC inlet
      3  compressor exit
      4  combustor exit before cooling air addition
      41 turbine  inlet after  cooling air addition
      45 HPT exit, LPT inlet
      49 LPT exit
      5  core nozzle
      6  core flow downstream
      7  fan nozzle
      8  fan flow downstream
"""
function tfsize!(gee, M0, T0, p0, a0, M2, M25,
      Feng, Phiinl, Kinl, iBLIc,
      BPR, pif, pilc, pihc,
      pid, pib, pifn, pitn,
      Ttf, ifuel, etab,
      epf0, eplc0, ephc0, epht0, eplt0,
      pifK, epfK,
      mofft, Pofft,
      Tt9, pt9, Tt4,
      epsl, epsh,
      icool,
      Mtexit, dTstrk, StA, efilm, tfilm,
      M4a, ruc,
      ncrowx, ncrow,
      epsrow, Tmrow,
      Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
      Δp_PreC, Δp_InterC, Δp_Regen)

      n = 6

      # from 'tfmap.inc'
      #        a     b     k     mo     da    c    d     C    D
      Cmapf = [3.50, 0.80, 0.03, 0.95, -0.50, 3.0, 6.0, 0.0, 0.0]
      Cmapl = [1.90, 1.00, 0.03, 0.95, -0.20, 3.0, 5.5, 0.0, 0.0]
      Cmaph = [1.75, 2.00, 0.03, 0.95, -0.35, 3.0, 5.0, 0.0, 0.0]

      #       Pcon   Ncon
      Tmapl = [0.15, 0.15]
      Tmaph = [0.15, 0.15]

      #---- fractional core mass flow convergence tolerance
      toler = 1.0e-12

      #---- mass offtake fraction update under-relaxation factor
      rlxfo = 0.8

      mcore = 0.0
      fo = 0.0
      Pom = 0.0

      #---- overall pressure ratio
      pic = pilc * pihc

      #
      # ===============================================================
      gas0 = IdealGases.Gas() #Initialize freestream gas

      #---- freestream static quantities
      gas0 = IdealGases.set_TP!(gas0, T0, p0)
      _, _, s0, dsdt, h0, dhdt, cp0, R0 = gas_unpack(gas0)
      
      u0 = M0 * a0

      # ===============================================================
      #---- freestream total quantities
      gast0 = deepcopy(gas0) #Initialize total freestream gas
      IdealGases.gas_Mach!(gast0, M0, 0.0)

      Tt0, pt0, st0, dsdt, ht0, dhdt, cpt0, Rt0 = gas_unpack(gast0)

      at0 = sqrt(Tt0 * Rt0 * cpt0 / (cpt0 - Rt0))

      # ===============================================================
      #---- offtake plume flow 9
      Trat = (p0 / pt9)^(Rt0 / cpt0)
      if (Trat < 1.0)
            u9 = sqrt(2.0 * cpt0 * Tt9 * (1.0 - Trat))
            rho9 = p0 / (Rt0 * Tt0 * Trat)
      else
            u9 = 0.0
            rho9 = p0 / (Rt0 * Tt0)
      end

      # ===============================================================
      #---- diffuser flow 0-2
      
      pt18 = pt0 * pid
      gast18 = deepcopy(gast0)
      gast18 = IdealGases.set_TP!(gast18, Tt0, pt18)

      Tt18, pt18, st18, _, ht18, _, cpt18, Rt18 = gas_unpack(gast0)

      #---- initial guesses for station 2 and 1.9
      pt2 = pt18
      Tt2 = Tt18
      pt19 = pt18
      Tt19 = Tt18


      if (Kinl == 0.0 && mofft == 0.0 && Pofft == 0.0)
            #----- single design pass will be sufficient
            npass = 1
      else
            #----- must use multiple passes to converge pt2,pt19 from inlet defect Kinl,
            #-     and to converge on offtake fractions fo, Pom
            npass = 60
      end

      sbfan = 0.0
      sbcore = 0.0
      for ipass = 1:npass

            # ===============================================================
            #---- set fan inlet conditions corrected for BLI
            if (ipass == 1)
                  #c      if(mcore == 0.0)
                  #----- don't know engine mass flow yet, so ignore any BLI mixing
                  if (iBLIc == 0)
                        sbfan = 0.0
                        sbcore = 0.0
                  else
                        sbfan = 0.0
                        sbcore = 0.0
                  end

            else
                  #----- account for inlet BLI defect via mass-averaged entropy
                  a2sq = at0^2 / (1.0 + 0.5 * (gam0 - 1.0) * M2^2)

                  if (iBLIc == 0)
                        #------ BL mixes with fan flow only
                        mmix = BPR * mcore * sqrt(Tt2 / Tt0) * pt0 / pt2
                        sbfan2 = Kinl * gam0 / (mmix * a2sq)
                        sbcore2 = 0.0
                  else
                        #------ BL mixes with fan + core flow
                        mmix = BPR * mcore * sqrt(Tt2 / Tt0) * pt0 / pt2 +
                               mcore * sqrt(Tt19 / Tt0) * pt0 / pt19
                        sbfan2 = Kinl * gam0 / (mmix * a2sq)
                        sbcore2 = sbfan
                  end

                  #----- update mixed-out entropies, with some underrelaxation       
                  rlxs = 0.85
                  sbfan = sbfan + rlxs * (sbfan2 - sbfan)
                  sbcore = sbcore + rlxs * (sbcore2 - sbcore)

            end

            #---- note: BL is assumed adiabatic, 
            #-     so Tt2,ht2,st2,cpt2,Rt2  will not change due to BL ingestion
            pt2 = pt18 * exp(-sbfan)
            gast2 = deepcopy(gast18)
            IdealGases.set_TP!(gast2, Tt18, pt2)
            Tt2, pt2, st2, _, ht2, _, cpt2, Rt2 = gas_unpack(gast2)

            pt19 = pt18 * exp(-sbcore)
            gast19 = deepcopy(gast18)
            IdealGases.set_TP!(gast19, Tt19, pt19)
            Tt19, pt19, st19, _, ht19, _, cpt19, Rt19 = gas_unpack(gast19)

            # ===============================================================
            #---- fan flow 2-7
            pifD = pif
            mbfD = 1.0
            mf = 1.0

            epf, epf_pf, epf_mf = ecmap(pif, mf, pifD, mbfD, Cmapf, epf0, pifK, epfK)

            gast21 = deepcopy(gast2)
            IdealGases.compress(gast21, pif, epf)

            Tt21, pt21, st21, _, ht21, _, cpt21, Rt21 = gas_unpack(gast21)

            #---- fan duct nozzle total quantities
            pt7 = pt21 * pifn
            
            gast7 = deepcopy(gast21)
            IdealGases.set_TP!(gast7, Tt21, pt7)
            Tt7, pt7, st7, _, ht7, _, cpt7, Rt7 = gas_unpack(gast7)

            # ===============================================================
            #---- Compressor precooler 19 - 19c
            pt19c = pt19 - Δp_PreC
            ht19c = ht19 + Δh_PreC
            
            gast19c = deepcopy(gast19)
            IdealGases.set_hP!(gast19c, ht19c, pt19c)
            Tt19c, pt19c, st19c, _, ht19c, _, cpt19c, Rt19c = gas_unpack(gast19c)

            #c      write(*,*) 'epf0 epf ', epf0, epf
            #
            # ===============================================================
            #---- LP compressor flow 19c - 25
            pilcD = pilc
            mblcD = 1.0
            ml = 1.0
            eplc, eplc_pl, eplc_ml = ecmap(pilc, ml, pilcD, mblcD, Cmapl, eplc0, 1.0, 0.0)

            gast25 = deepcopy(gast19c)
            IdealGases.compress(gast25, pilc, eplc0)

            Tt25, pt25, st25, _, ht25, _, cpt25, Rt25 = gas_unpack(gast25)

            # ===============================================================
            #---- Compressor intercooler 25 - 25c
            pt25c = pt25 - Δp_InterC
            ht25c = ht25 + Δh_InterC
            
            gast25c = deepcopy(gast25)
            IdealGases.set_hP!(gast25c, ht25c, pt25c)
            Tt25c, pt25c, st25c, _, ht25c, _, cpt25c, Rt25c = gas_unpack(gast25c)

            # ===============================================================
            #---- HP compressor flow 25c - 3
            pihcD = pihc
            mbhcD = 1.0
            mh = 1.0
            ephc, ephc_ph, ephc_mh = ecmap(pihc, mh, pihcD, mbhcD, Cmaph, ephc0, 1.0, 0.0)

            gast3 = deepcopy(gast25c)
            IdealGases.compress(gast3, pihc, ephc0)

            Tt3, pt3, st3, _, ht3, _, cpt3, Rt3 = gas_unpack(gast3)


            # ===============================================================
            #---- combustor flow 3-4   (ffb = mdot_fuel/mdot_burner)
            fuel = "Jet-A(g)"
            ffb, gast4 = IdealGases.gas_burn(gast3, fuel, Ttf, Tt4, etab, 0.0)
            
            pt4 = pt3 * pib
            IdealGases.set_TP!(gast4, Tt4, pt4)

            _, _, st4, _, ht4, _, cpt4, Rt4 = gas_unpack(gast4)

            # ===============================================================

            if (icool == 0)
                  #----- no cooling air present... station 41 is same as 4
                  ff = ffb * (1.0 - fo)

                  gast41 = deepcopy(gast4)
                  Tt41, pt41, st41, _, ht41, _, cpt41, Rt41 = gas_unpack(gast41)

                  #
                  #----------------------------------------------------------------
            else
                  #----- cooling air is present... calculate station 41

                  #----- hot-section temperature ratio for each blade row (for cooling model)
                  gmi4 = Rt4 / (cpt4 - Rt4)
                  Trrat = 1.0 / (1.0 + 0.5 * gmi4 * Mtexit^2)

                  # Heat exchanger to cool turbine cooling air
                  ht_tc = ht3 + Δh_TurbC #Specific enthalpy of turbine cooling air
                  gast_tc = deepcopy(gast3)
                  IdealGases.set_hP!(gast_tc, ht_tc, gast3.P)
                  Tt_tc = gast_tc.T

                  if (icool == 1)
                        #------ epsrow(.) is assumed to be passed in.. calculate Tmrow(.)
                        Tmrow = Tmcalc(ncrowx, ncrow,
                        Tt_tc, Tt4, dTstrk, Trrat,
                              efilm, tfilm, StA, epsrow)
                  else
                        #------ calculate cooling mass flow ratios epsrow(.) to get specified Tmrow(.)
                        ncrow, epsrow, epsrow_Tt3, epsrow_Tt4, epsrow_Trr = mcool(ncrowx,
                              Tmrow, Tt_tc, Tt4, dTstrk, Trrat,
                              efilm, tfilm, StA)
                  end

                  #----- total cooling-flow fraction
                  fc = 0.0
                  for icrow = 1:ncrow
                        fc = fc + (1.0 - fo) * epsrow[icrow]
                  end

                  if (fc >= 0.99)
                        error("TFSIZE: Excessive cooling flow", 
                              "\n\tmcool/mcore = ", fc, 
                              "\n\tTt3 Tt4 Tmetal ", Tt3, " K, ", Tt4, " K, ", Tmrow[1], " K")
                  end

                  #----- set ff = mdot_fuel/mdot_core = ffb * mdot_burner/mdot_core
                  ff = (1.0 - fo - fc) * ffb


                  #----- speed at start-of-mixing station 4a
                  gas4a = deepcopy(gast4)
                  IdealGases.gas_Mach!(gas4a, 0.0, M4a)

                  T4a, p4a, s4a, _, h4a, _, cp4a, R4a = gas_unpack(gas4a)

                  u4sq = max(2.0 * (ht4 - h4a), 0.0)
                  u4a = sqrt(u4sq)

                  #----- exit speed of cooling air at station 4a
                  uc = ruc * u4a

                  #----- IGV exit mixing
                  frac4 = (1.0 - fo - fc + ff) / (1.0 - fo + ff)
                  fracm = fc / (1.0 - fo + ff)

                  #----- mixed constituent fraction vector from mass equation
                  # for i = 1:nair
                  #       lambdap[i] = frac4 * lambda[i] + fracm * alpha[i]
                  # end

                  gas41 = IdealGases.gas_mixing(gast3, gast4, frac4 / fracm)

                  #----- mixed total enthalpy from enthalpy equation
                  ht41 = frac4 * ht4 + fracm * ht3

                  #----- mixed velocity from momentum equation, assuming constant static pressure
                  p41 = p4a
                  u41 = frac4 * u4a + fracm * uc

                  #----- static temperature from static enthalpy
                  h41 = ht41 - 0.5 * u41^2

                  IdealGases.set_hP!(gas41, h41, p41)

                  #----- all static quantities, from static temperature
                  T41, _, s41, _, h41, _, cp41, R41 = gas_unpack(gas41)

                  #----- all stagnation quantities, from total-static enthalpy difference
                  dhb = ht41 - h41
                  epi = 1.0
                  gast41 = IdealGases.gas_Deltah(gas41, dhb, epi)

                  Tt41, pt41, st41, _, ht41, _, cpt41, Rt41 = gas_unpack(gast41)

            end

            # ===============================================================
            #---- LPT and HPT work, per unit mass flow
            dhfac = -(1.0 - fo) / (1.0 - fo + ff) / (1.0 - epsh)
            dlfac = -1.0 / (1.0 - fo + ff) / (1.0 - epsl)

            dhht = (ht3 - ht25c) * dhfac
            dhlt = (ht25 - ht19c + BPR * (ht21 - ht2) + Pom) * dlfac

            #---- HPT flow
            #     Trh =  Tt41/(Tt41 + dhht/cpt41)
            #     gexh = cpt41/(Rt41*epht0)
            #     pihtD = Trh^gexh

            epi = 1.0 / epht

            gast45 = IdealGases.gas_Deltah(gast41, dhht, epi)
            Tt45, pt45, st45, _, ht45, _, cpt45, Rt45 = gas_unpack(gast45)


            epi = 1.0 / eplt
            gast49 = IdealGases.gas_Deltah(gast45, dhlt, epi)
            Tt49, pt49, st49, _, ht49, _, cpt49, Rt49 = gas_unpack(gast49)

            # ===============================================================
            #---- Regenerative cooling heat exchanger 49 - 49c
            pt49c = pt49 - Δp_Regen
            ht49c = ht49 + Δh_Regen

            gast49c = deepcopy(gast49)
            IdealGases.set_hP!(gast49c, ht49c, pt49c)
            Tt49c, pt49c, st49c, _, ht49c, _, cpt49c, Rt49c = gas_unpack(gast49c)

            # ===============================================================
            #---- Turbine nozzle 49c - 5

            gast5 = deepcopy(gast49)
            pt5 = pt49c * pitn
            IdealGases.set_TP!(gast5, Tt49c, pt5)
            Tt5, pt5, st5, _, ht5, _, cpt5, Rt5 = gas_unpack(gast5)

            #
            # ===============================================================
            #---- fan plume flow 7-8, use alpha mass fraction (air)
            gast8 = deepcopy(gast7)
            Tt8, pt8, st8, _, ht8, _, cpt8, Rt8 = gas_unpack(gast8)

            pratfn = p0 / pt8

            gas8 = IdealGases.expand(gast8, pratfn, 1.0)
            T8, p8, s8, _, h8, _, cp8, R8= gas_unpack(gas8)
            if (h8 >= ht8)

                  Lconv = false

                  u8 = 0.001 * sqrt(R8 * T8)
                  error("TFSIZE: Negative fan plume velocity", 
                        "\n\tpt2,  Tt2  = ", pt2, " Pa, ",  Tt2, " K",
                        "\n\tpt8,  Tt8  = ", pt8, " Pa, ",  Tt8, " K", 
                        "\n\tp8,  T8  = "  , p8,  " Pa, ",  T8,  " K", 
                        "\n\tpif,  BPR  = ", pif, " Pa, ",  BPR)
            else
                  u8 = sqrt(2.0 * (ht8 - h8))
            end
            rho8 = gas8.rho

            # ===============================================================
            #---- core plume flow 5-6, use lambdap mass fraction (combustion products)
            gast6 = deepcopy(gast5)
            Tt6, pt6, st6, _, ht6, _, cpt6, Rt6 = gas_unpack(gast6)

            prattn = p0 / pt6
            gas6 = IdealGases.expand(gast6, prattn, 1.0)
            T6, p6, s6, _, h6, _, cp6, R6 = gas_unpack(gas6)
            
            if (h6 >= ht6)
                  
                  Lconv = false
                  u6 = 0.001 * sqrt(R6 * T6)
                  error("TFSIZE: Negative core plume velocity", 
                        "\n\tpt2,  Tt2  = ", pt2, " Pa, ",  Tt2, " K",
                        "\n\tpt3,  Tt3  = ", pt3, " Pa, ",  Tt3, " K", 
                        "\n\tpt4,  Tt4  = ", pt4, " Pa, ",  Tt4, " K", 
                      "\n\tpt41,  Tt41  = ", pt41," Pa, ",  Tt41," K", 
                        "\n\tpt6,  Tt6  = ", pt6, " Pa, ",  Tt6, " K", 
                        "\n\tp6,  T6  = "  , p6,  " Pa, ",  T6,  " K", 
                        "\n\tpif,  BPR  = ", pif, " Pa, ",  BPR)
            else
                  u6 = sqrt(2.0 * (ht6 - h6))
            end

            rho6 = gas6.rho

            #      write(*,*) 'Pt6   u6 ', pt6, u6
            #
            # ===============================================================
            #---- effective fuel heating value, over states 3, 4  (just for info)
            cpa = 0.5 * (cpt3 + cpt4)
            hfuel = cpa * (Tt4 - Tt3 + ffb * (Tt4 - Ttf)) / (etab * ffb)

            #---- effective fuel heating value, over states 3, 4.1  (just for info)
            #      cpa = 0.5*(cpt3+cpt41)
            #      hfuel = cpa*((1.0-fo)*(Tt41-Tt3) + ff*(Tt41-Ttf)) / (etab*ff)

            # ===============================================================
            #---- size core mass flow

            #---- store current values for better update, convergence checks
            mcold = mcore
            foold = fo

            #---- added effective net thrust from dissipation in ingested streamtube
            if (u0 == 0.0)
                  Finl = 0.0
            else
                  Finl = Phiinl / u0
            end

            #---- set core mass flow from specified effective net thrust
            mcore = (Feng - Finl) /
                    ((1.0 - fo + ff) * u6 - u0 + BPR * (u8 - u0) + fo * u9)

            #---- corresponding new offtake mass flow fraction
            fonew = mofft / mcore
            dfo = fonew - foold

            fo = fo + rlxfo * dfo

            #---- estimate better new mass flow, compensating for change in mass offtake
            mfac = min(2.0, 1.0 / (1.0 - dfo))
            mcore = mcore * mfac

            #---- power offtake per mass flow
            Pom = Pofft / mcore

            #---- overall Fsp and TSFC
            Fsp = Feng / (u0 * mcore * (1.0 + BPR))
            if (Feng <= 0.0)
                  TSFC = 0.0
            else
                  TSFC = (gee * ff * mcore) / Feng
            end

            # ===============================================================

            M8 = u8 / sqrt(cp8 * R8 / (cp8 - R8) * T8)
            if (M8 < 1.0)
                  #----- subsonic fan plume... fan nozzle flow is same as plume
                  gas7 = deepcopy(gas8)
                  T7, p7, s7, _, h7, _, cp7, R7 = gas_unpack(gas7)
                  u7 = u8
            else
                  #----- supersonic fan plume... fan nozzle is choked
                  M7 = 1.0

                  gas7 = deepcopy(gast7)
                  IdealGases.gas_Mach!(gas7, 0.0, M7, 1.0)
                  T7, p7, s7, _, h7, _, cp7, R7 = gas_unpack(gas7)
                  u7 = sqrt(2.0 * (ht7 - h7))
            end
            rho7 = gas7.rho

            #---- size fan  nozzle and plume areas
            A7 = BPR * mcore / (rho7 * u7)
            A8 = BPR * mcore / (rho8 * u8)

            # ===============================================================
            M6 = u6 / sqrt(cp6 * R6 / (cp6 - R6) * T6)
            #      write(*,*) 'u6,M6', u6,M6
            if (M6 < 1.0)
                  #----- subsonic core plume... core nozzle flow is same as plume
                  gas5 = deepcopy(gas6)
                  T5, p5, s5, _, h5, _, cp5, R5 = gas_unpack(gas5)
                  u5 = u6
                  
            else
                  #----- supersonic core plume... core nozzle is choked
                  M5 = 1.0

                  gas5 = deepcopy(gast6)
                  IdealGases.gas_Mach!(gas5, 0.0, M5, 1.0)
                  T5, p5, s5, _, h5, _, cp5, R5 = gas_unpack(gas5)
                  u5 = sqrt(2.0 * (ht5 - h5))
                  
            end

            rho5 = gas5.rho

            #---- size core nozzle and plume areas
            A5 = (1.0 - fo + ff) * mcore / (rho5 * u5)
            A6 = (1.0 - fo + ff) * mcore / (rho6 * u6)

            if (u9 == 0.0)
                  A9 = 0.0
            else
                  A9 = fo * mcore / (rho9 * u9)
            end

            # ===============================================================
            #---- size fan and compressor areas
            gas2 = deepcopy(gast2)
            IdealGases.gas_Mach!(gas2, 0.0, M2, 1.0)
            T2, p2, s2, _, h2, _, cp2, R2 = gas_unpack(gas2)
            rho2 = gas2.rho
            u2 = sqrt(2.0 * (ht2 - h2))

            gas19c = deepcopy(gast19c)
            IdealGases.gas_Mach!(gas19c, 0.0, M2, 1.0)
            T19c, p19c, s19c, _, h19c, _, cp19c, R19c = gas_unpack(gas19c)
            rho19c = gas19c.rho
            u19c = sqrt(2.0 * (ht19c - h19c))

            A2 = BPR * mcore / (rho2 * u2) + mcore / (rho19c * u19c)

            gas25c = deepcopy(gast25c)
            IdealGases.gas_Mach!(gas25c, 0.0, M25, 1.0)
            T25c, p25c, s25c, _, h25c, _, cp25c, R25c = gas_unpack(gas25c)
            rho25c = gas25c.rho
            u25c = sqrt(2.0 * (ht25c - h25c))

            A25 = (1.0 - fo) * mcore / (rho25c * u25c)

            if (ipass >= 2)
                  dmfrac = 1.0 - mcold / mcore

                  if (abs(dmfrac) < toler)

                        # ===============================================================
                        #---- calculate component efficiencies  (informative only -- not needed here)
                        etaf = 0.0
                        etalc = 0.0
                        etahc = 0.0
                        etaht = 0.0
                        etalt = 0.0

                        #---- fan
                        gast21i = IdealGases.PressureRatio(gast2, pif, 1.0)
                        Tt21i, pt21i, st21i, _, ht21i, _, cpt21i, Rt21i = gas_unpack(gast21i)
                        etaf = (ht21i - ht2) / (ht21 - ht2)

                        #---- LP compressor
                        gast25i = IdealGases.PressureRatio(gast19c, pilc, 1.0)
                        Tt25i, pt25i, st25i, _, ht25i, _, cpt25i, Rt25i = gas_unpack(gast25i)
                        etalc = (ht25i - ht19c) / (ht25 - ht19c)

                        #---- HP compressor
                        gast3i = IdealGases.PressureRatio(gast25c, pihc, 1.0)
                        Tt3i, pt3i, st3i, _, ht3i, _, cpt3i, Rt3i = gas_unpack(gast3i)
                        etahc = (ht3i - ht25c) / (ht3 - ht25c)

                        #---- HP turbine
                        piht = pt45 / pt41

                        gast45i = IdealGases.PressureRatio(gast41, piht, 1.0)
                        Tt45i, pt45i, st45i, _, ht45i, _, cpt45i, Rt45i = gas_unpack(gast45i)
                        etaht = (ht45 - ht41) / (ht45i - ht41)

                        #---- LP turbine
                        pilt = pt49 / pt45
                        gast49i = IdealGases.PressureRatio(gast45, pilt, 1.0)
                        Tt49i, pt49i, st49i, _, ht49i, _, cpt49i, Rt49i = gas_unpack(gast49i)
                        etalt = (ht49 - ht45) / (ht49i - ht45)
                        
                        Lconv = true
                        return epsrow, Tmrow,
                        TSFC, Fsp, hfuel, ff, mcore,
                        Tt0, ht0, pt0, cpt0, Rt0,
                        Tt18, ht18, pt18, cpt18, Rt18,
                        Tt19, ht19, pt19, cpt19, Rt19,
                        Tt19c, ht19c, pt19c, cpt19c, Rt19c,
                        Tt2, ht2, pt2, cpt2, Rt2,
                        Tt21, ht21, pt21, cpt21, Rt21,
                        Tt25, ht25, pt25, cpt25, Rt25,
                        Tt25c, ht25c, pt25c, cpt25c, Rt25c,
                        Tt3, ht3, pt3, cpt3, Rt3,
                        ht4, pt4, cpt4, Rt4,
                        Tt41, ht41, pt41, cpt41, Rt41,
                        Tt45, ht45, pt45, cpt45, Rt45,
                        Tt49, ht49, pt49, cpt49, Rt49,
                        Tt5, ht5, pt5, cpt5, Rt5,
                        Tt7, ht7, pt7, cpt7, Rt7,
                        u0,
                        T2, u2, p2, cp2, R2, A2,
                        T25c, u25c, p25c, cp25c, R25c, A25,
                        T5, u5, p5, cp5, R5, A5,
                        T6, u6, p6, cp6, R6, A6,
                        T7, u7, p7, cp7, R7, A7,
                        T8, u8, p8, cp8, R8, A8,
                        u9, A9,
                        epf, eplc, ephc, epht, eplt,
                        etaf, etalc, etahc, etaht, etalt,
                        Lconv
                  end

            end


      end

      if (npass > 1)
            println("TFSIZE: Convergence failed.  dm/m = ", dmfrac)
      end

end # tfsize