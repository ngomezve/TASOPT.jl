@testset "structural components" begin

  # Fuselage
  include(joinpath(TASOPT.__TASOPTroot__, "../test/default_structures.jl"))
  fuselage = ac_test.fuselage

  Nland, Wpaymax, Wengtail, 
                          nftanks,
                          Waftfuel,  Wftank_single, ltank, xftank_fuse, tank_placement,
                          Δp,
                          Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
                          bv, λv, nvtail,
                          xhtail, xvtail,
                          xwing, xwbox, cbox,
                          xeng = [6.0, 219964.5779, 0.0, 0, 0.0, 0.0, 0.0, 0.0, "", 56016.1756316446, 14400.815502104786, 
        9645.22120844135, 0.4, 0.7, 1.2379754637499652e6, 1.0469704911343623e6, 7.958844386073947, 
        0.3, 1.0, 36.212475021199154, 35.05057234350367, 18.465770971525686, 15.688410884477133, 
        3.119635043953234, 14.164410884477133 ]
  
  cabVol = TASOPT.fusew!(fuselage, Nland, Wpaymax, Wengtail, 
                          nftanks,
                          Waftfuel,  Wftank_single, ltank, xftank_fuse, tank_placement,
                          Δp,
                          Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
                          bv, λv, nvtail,
                          xhtail, xvtail,
                          xwing, xwbox, cbox,
                          xeng)
    
    out = [fuselage.skin.thickness, fuselage.cone.thickness, 
            fuselage.layout.thickness_webs, fuselage.floor.thickness, 
            fuselage.bendingmaterial_h.weight.x, fuselage.bendingmaterial_v.weight.x,
            fuselage.shell.EIh, fuselage.bendingmaterial_h.EIh,
            fuselage.shell.EIv, fuselage.bendingmaterial_v.EIh,
            fuselage.shell.GJ , fuselage.cone.GJ,
            fuselage.shell.weight.W, fuselage.cone.weight.W, 
            fuselage.window.W, fuselage.insulation.W, 
            fuselage.floor.weight.W,
            fuselage.bendingmaterial_h.weight.W, fuselage.bendingmaterial_v.weight.W,
            fuselage.weight,fuselage.moment,
            cabVol]

  test_out = [0.0010592916489634418, 0.0012236016395195024, 0.0, 0.0009061565981319232, 
  31.85585727200865, 31.098482101397437, 2.361320105799576e9, 2.4685708558261826e10, 
  2.317349483950291e9, 2.4685708558261826e10, 1.5708404995400493e9, 1.6421853018493836e9, 
  20947.352673154768, 5667.810474837753, 11269.980000000001, 5018.7283576257105, 13415.796495703155, 
  14314.082074792976, 7380.925148634968, 198042.16150624937, 3.620261283824471e6, 414.54271063740043]
    
@test all(isapprox.(out, test_out))

# end Fuselage weight

# calculate_centroid_offset:
  b = 35.486921629195265
  bs = 10.113772664320649
  bo = 3.6067999999999998
  lambdat = 0.25000000000000000
  lambdas = 0.69999999999999996
  sweep = 26.000000000000000
  dxwing = 2.5792921907840762
  macco = 0.68475212288241949

dx_out, macco_out = TASOPT.structures.calculate_centroid_offset(b,bs,bo,lambdat,lambdas,sweep)
@test dx_out ≈ dxwing
@test macco_out ≈ macco 

# end calculate_centroid_offset


#surfw:
  po = 114119.45308868506
  gammat = 0.225
  gammas = 0.8665999999999999
  Nlift = 3.0000000000000000
  Weng1 = 30083.769249783934
  fLt = -0.05
  sigfac = 1.0
  rhofuel = 817.00000000000000

  fort_Wwing = 129055.8385739653
  fort_Wsinn = 11643.760847427178
  fort_Wsout = 20444.306372488714
  fort_dyWsinn = 18266.601840430056
  fort_dyWsout = 96970.09656550126
  fort_Wfcen =  30488.07253761438
  fort_Wfinn = 38682.37013612792
  fort_Wfout = 46796.44312114958
  fort_dxWfinn = 29597.79354199928 
  fort_dxWfout = 189307.87250917125
  fort_dyWfinn =  60684.46980143201
  fort_dyWfout = 221961.82769430208
  fort_lstrutp = 0.0


Wwing,Wsinn,Wsout,
        dyWsinn,dyWsout,
        Wfcen,Wfinn,Wfout,
        dxWfinn,dxWfout,
        dyWfinn,dyWfout,lstrutp = TASOPT.surfw!(wing, po, gammat, gammas,
                                            Nlift, Weng1, 0, 0.0, 0, 0.0,
                                            sigfac, rhofuel)


@test fort_Wwing ≈ Wwing 
@test fort_Wsinn ≈ Wsinn 
@test fort_Wsout ≈ Wsout 
@test fort_dyWsinn ≈ dyWsinn 
@test fort_dyWsout ≈ dyWsout 
@test fort_Wfcen ≈ Wfcen 
@test fort_Wfinn ≈ Wfinn 
@test fort_Wfout ≈ Wfout 
@test fort_dxWfinn ≈ dxWfinn 
@test fort_dxWfout ≈ dxWfout 
@test fort_dyWfinn ≈ dyWfinn 
@test fort_dyWfout ≈ dyWfout 
@test fort_lstrutp ≈ lstrutp 
#end surfw

#surfw for Htail
poh = 115893.98734144184
λhs = 1.0
fLt = -0.05
tauwebh = 1.378913257881327e8
σcaph = 2.0684848484848484e8
surft_f_out = [14400.81547163942, 14069.611170000926, 0.0011568849664072272, 0.0023905578555627194, 1.896322960387795e8, 1.2616774558497725e9, 1.794003547786858e8]

TASOPT.surfw!(htail, poh, htail.outboard.λ, λhs,
0.0, 0.0, 0, 0.0, 0, 0.0,
sigfac, rhofuel)

surft_out = [htail.weight, htail.dxW, htail.outboard.webs.thickness, htail.outboard.caps.thickness, htail.outboard.EI[1], htail.outboard.EI[4], htail.outboard.GJ]

@test all(isapprox.(surft_out, surft_f_out))
#end surft


#tailpo:
  Sh = 42.443587259163984
  ARh = 6.0000000000000000
  lambdah = 0.25000000000000000
  qne = 12692.519555311032
  CLhmax = 2.0000000000000000
  fort_bh = 15.958117794871169
  fort_coh = 4.2554980786323124
  fort_poh = 108025.98516125829

poh = TASOPT.structures.tailpo!(htail, Sh, qne)

@test fort_bh ≈ htail.layout.span
@test fort_coh ≈ htail.layout.root_chord
@test fort_poh ≈ poh
#end tailpo:

# #Test fuselage layout updating for tank

parg = zeros(igtotal)
fuselage.layout.x_start_cylinder = 6.096
fuselage.layout.x_end_cylinder = 29.5656
fuselage.layout.x_pressure_shell_aft = 31.0896
parg[igdxeng2wbox] = 1.5239999999999991
fuselage.APU.r = [36.576, 0.0, 0.0]
fuselage.layout.x_end = 37.7952
fuselage.layout.x_cone_end = 35.6616
htail.layout.box_x = 34.8996
vtail.layout.box_x = 33.528
wing.layout.box_x = 16.04432532088372
parg[igxeng] = wing.layout.box_x - parg[igdxeng2wbox]
fuselage.layout.x_cone_end = fuselage.layout.x_cone_end * 0.52484 

pari = zeros(iitotal)
pari[iinftanks] = 1

#Update fuel tank length and check changes
parg[iglftank] = 5.0
TASOPT.update_fuse!(fuselage, wing, htail, vtail, pari, parg)

update_fuse_out = [fuselage.layout.x_end_cylinder, 
fuselage.layout.x_pressure_shell_aft, 
fuselage.layout.x_cone_end, 
fuselage.APU.x,
fuselage.layout.x_end,
fuselage.HPE_sys.x,
htail.layout.box_x,
vtail.layout.box_x,
parg[igxeng]]

update_fuse_out_test = [35.175200000000004, 36.699200000000005, 24.326234144000004, 42.18560000000001, 43.40480000000001, 12.767380728136962, 40.50920000000001, 39.137600000000006, 14.52032532088372]
@test all(isapprox.(update_fuse_out, update_fuse_out_test))


#Return to original points?
pari[iinftanks] = 0.0
parg[iglftank] = 0.0
TASOPT.update_fuse!(fuselage, wing, htail, vtail, pari, parg)

update_fuse_out = [fuselage.layout.x_end_cylinder, 
fuselage.layout.x_pressure_shell_aft, 
fuselage.layout.x_cone_end, 
fuselage.APU.x,
fuselage.layout.x_end,
fuselage.HPE_sys.x,
htail.layout.box_x,
vtail.layout.box_x,
parg[igxeng]]

update_fuse_out_test = [29.5656, 31.0896, 18.716634144, 36.57600000000001, 37.79520000000001, 9.82323826413696, 34.89960000000001, 33.528000000000006, 14.52032532088372]
@test all(isapprox.(update_fuse_out, update_fuse_out_test))

#Test cabin resizing
parg = zeros(igtotal)
fuselage.layout.x_start_cylinder = 6.096
fuselage.layout.x_end_cylinder = 29.5656
fuselage.layout.x_pressure_shell_aft = 31.0896
fuselage.layout.l_cabin_cylinder = 23.4696
parg[igdxeng2wbox] = 1.5239999999999991
fuselage.APU.r = [36.576, 0.0, 0.0]
fuselage.layout.x_end = 37.7952
fuselage.layout.x_cone_end = 35.6616
parg[igxhbox ] = 34.8996
parg[igxvbox ] = 33.528
parg[igxwbox] = 16.04432532088372
parg[igxeng] = parg[igxwbox] - parg[igdxeng2wbox]
fuselage.layout.x_cone_end = fuselage.layout.x_cone_end * 0.52484 

parg[igseatpitch] = 0.762
parg[igseatwidth] = 0.4826
parg[igaislehalfwidth] = 0.254
parg[igWpaymax] = 219964.5779
fuselage.layout.cross_section.radius = 2.5 #Change radius to 2.5 m

pari = zeros(iitotal)
pari[iidoubledeck] = 0

fuse_tank = TASOPT.fuselage_tank()

TASOPT.update_fuse_for_pax!(pari, parg, fuselage, fuse_tank)

parg_check = [14.584924835954398, 16.108924835954397, 35.05200000000001, 33.680400000000006, 219964.5779, 1.5239999999999991, 0.762, 0.4826, 0.254]
parg_nz = deepcopy(parg)
deleteat!(parg_nz, parg_nz .== 0)
for (i,item) in enumerate(parg_nz) #For every nonzero element in parg
  @test parg_nz[i] ≈ parg_check[i]
end

@test fuselage.layout.x_pressure_shell_aft ≈ 31.24200000000001
@test fuselage.layout.x_cone_end ≈ 18.86903414400001
@test fuselage.layout.x_end ≈ 37.94760000000001

#Test minimum radius calculation
Rmin = TASOPT.structures.find_minimum_radius_for_seats_per_row(5, ac_test)
@test Rmin ≈ 1.7113052179793784
end