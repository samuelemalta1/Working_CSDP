import os
import sys
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtWidgets import QDialog, QApplication, QWidget, QMainWindow, QLabel, QDialogButtonBox
from PyQt5.QtWidgets import QMessageBox
import main
import Aux_classes as aux

set_images_path = "../Turbomachinery_code/"
class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        loadUi("GUI.ui",self)

        # Initialize Font
        screen = QApplication.primaryScreen()
        size = screen.size()
        fontSize = int(size.width()/100)
        self.setStyleSheet(''' font-size: ''' + str(fontSize) + '''px; ''')

        # Initialize graphics


        #tabs
        self.tab_inputs.setCurrentIndex(0)
        self.tab_outputs.setCurrentIndex(0)

        # Connects
        #Main
        self.line_thrust.editingFinished.connect(self.checkThrust);  self.checkThrust()
        self.line_time.editingFinished.connect(self.checktime); self.checktime()
        self.line_pa.editingFinished.connect(self.checkPa); self.checkPa()
        self.line_O_F.editingFinished.connect(self.checkOF); self.checkOF()
        self.line_Otank_pres.editingFinished.connect(self.checkOtankPres); self.checkOtankPres()
        self.line_Ftank_pres.editingFinished.connect(self.checkFtankPres); self.checkFtankPres()

        #Nozzle
        self.nozzle_changed(0)
        self.line_De_max.editingFinished.connect(self.checkDeMax);  self.checkDeMax()
        self.line_eps_max.editingFinished.connect(self.checkEpsMax);  self.checkEpsMax()
        self.line_conv_ang.editingFinished.connect(self.checkConvAng);  self.checkConvAng()
        self.line_div_ang.editingFinished.connect(self.checkDivAng);  self.checkDivAng()
        self.line_par_0ang.editingFinished.connect(self.checkParAng0);  self.checkParAng0()
        self.line_par_fang.editingFinished.connect(self.checkParAngf);  self.checkParAngf()
        self.line_curv_conv.editingFinished.connect(self.checkCurvConv);  self.checkCurvConv()
        self.line_curv_div.editingFinished.connect(self.checkCurvDiv);  self.checkCurvDiv()

        #Combustion chamber
        self.line_SF.editingFinished.connect(self.checkSF);  self.checkSF()

        #Turbo
        self.cycle_changed(0)
        self.line_Opump_eff.editingFinished.connect(self.checkOpumpEff);  self.checkOpumpEff()
        self.line_Fpump_eff.editingFinished.connect(self.checkFpumpEff);  self.checkFpumpEff()
        self.line_turb_eff.editingFinished.connect(self.checkTurbEff);  self.checkTurbEff()
        self.line_mech_eff.editingFinished.connect(self.checkMechEff);  self.checkMechEff()
        self.line_Ele_power.editingFinished.connect(self.checkElePow);  self.checkElePow()
        self.line_valve_losses.editingFinished.connect(self.checkValveLoss);  self.checkValveLoss()

        #Cooling


        #Injectors
        self.line_Cd.editingFinished.connect(self.checkCd);  self.checkCd()

        #Ignitors
        self.line_tign.editingFinished.connect(self.checktIgn);  self.checktIgn()

        #Materials


    #Run
    def Run(self):
        #get propellant, cycle, nozzle, materials, etc...
        main.prop = aux.Propellant(self.combo_prop.currentIndex())
        main.default.cycle_type = self.combo_cycle.currentIndex()

        #run main
        check = main.Main(main.dat)
        if(not check): 
            msg = QMessageBox()
            msg.setWindowTitle("Calculation error!")
            msg.setText("A fatal error occured and the calculation was aborted")
            msg.exec_()

        #output data
        self.line_Isp.setText(str(main.dat.Isp))
        self.line_cstar.setText(str(main.dat.cstar))
        self.line_m.setText(str(main.dat.m))
        self.line_prop_mass.setText(str(main.dat.Mprop))
        self.line_mass.setText(str(main.dat.Mtot))
        self.line_cost.setText(str(main.dat.cost))
        self.line_reliability.setText(str(main.dat.rel))

        #Nozzle
        self.line_Dt.setText(str(main.dat.Dt))
        self.line_De.setText(str(main.dat.De))
        self.line_eps.setText(str(main.dat.Eps))
        self.line_len.setText(str(main.dat.L_total))
        self.line_len_conv.setText(str(main.dat.L_con))
        self.line_len_div.setText(str(main.dat.L_div))
        self.line_m_nozz.setText(str(main.dat.m_nozz))
        self.line_nozz_mass.setText(str(main.dat.Mnoz))

        #Combustion chamber
        self.line_Pc.setText(str(main.dat.Pc))
        self.line_Tc.setText(str(main.dat.Tc))
        self.line_Dc.setText(str(main.dat.Dc))
        self.line_Lc.setText(str(main.dat.Chamber_L))
        self.line_tc.setText(str(main.dat.ThicknessChamber))
        self.line_Mc.setText(str(main.dat.chambermass))

        #turbo
        self.line_Opump_power.setText(str(main.dat.W_Opump))
        self.line_Fpump_power.setText(str(main.dat.W_Fpump))
        self.line_Turb_power.setText(str(main.dat.W_turb))
        self.line_P_inj_in.setText(str(main.dat.ptinj))

        #Cooling


        #Injectors
        self.line_n_oxinj.setText(str(main.dat.n_ox))
        self.line_n_finj.setText(str(main.dat.n_f))
        self.line_vel_ox.setText(str(main.dat.v_iox))
        self.line_vel_fuel.setText(str(main.dat.v_if))

        #Ignitors
        self.line_Mign.setText(str(main.dat.Igniter_compound))


        #Materials


    #Blocks
    def OptOF(self, ischeck : bool):
        if ischeck:
            main.dat_O_F = 0
            self.line_O_F.setEnabled(False)
        else:
            self.line_O_F.setEnabled(True)

    def cycle_changed(self,i : int):
        if i == 4:
            self.line_Ele_power.setEnabled(True)
            self.line_turb_eff.setEnabled(False)
        else:
            self.line_Ele_power.setEnabled(False)
            if i == 5:
                self.line_Opump_eff.setEnabled(False)
                self.line_Fpump_eff.setEnabled(False)
                self.line_turb_eff.setEnabled(False)
                self.line_mech_eff.setEnabled(False)
            else:
                self.line_Opump_eff.setEnabled(True)
                self.line_Fpump_eff.setEnabled(True)
                self.line_turb_eff.setEnabled(True)
                self.line_mech_eff.setEnabled(True)

    def nozzle_changed(self,i : int):
        if i==0:
            main.default.Nozzle_type = 0
            self.line_par_0ang.setEnabled(False)
            self.line_par_fang.setEnabled(False)
            self.line_div_ang.setEnabled(True)
        else:
            main.default.Nozzle_type = 1
            self.line_par_0ang.setEnabled(True)
            self.line_par_fang.setEnabled(True)
            self.line_div_ang.setEnabled(False)

    #Checks
    def checkThrust(self):
        thrust = float(self.line_thrust.text())
        if thrust > 0.0 and thrust < 1.0e10:
            main.dat.Thrust = thrust;
        else:
            self.line_thrust.setText("15000.0")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid thrust, try again.")
            msg.exec_()

    def checktime(self):
        var = float(self.line_time.text())
        if var > 0.0 and var < 1000:
            main.dat.time = var;
        else:
            self.line_time.setText("10.0")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid thrust time, try again.")
            msg.exec_()

    def checkPa(self):
        var = float(self.line_pa.text())
        if var > 0.0 and var < 1.0e6:
            main.dat.Pa = var;
        else:
            self.line_pa.setText("101325.0")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid ambient pressure, try again.")
            msg.exec_()

    def checkOF(self):
        var = float(self.line_O_F.text())
        if var > 0.0 and var < 100:
            main.dat.O_F = var;
        else:
            self.line_O_F.setText("5.0")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid O/F, try again.")
            msg.exec_()

    def checkOtankPres(self):
        var = float(self.line_Otank_pres.text())
        if var > 0.0 and var < 1.0e10:
            main.default.p_to = var;
        else:
            self.line_Otank_pres.setText("1.0e5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid oxidizer tank pressure, try again.")
            msg.exec_()

    def checkFtankPres(self):
        var = float(self.line_Ftank_pres.text())
        if var > 0.0 and var < 1.0e10:
            main.default.ptf = var;
        else:
            self.line_Ftank_pres.setText("1.0e5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid fuel tank pressure, try again.")
            msg.exec_()

    def checkOpumpEff(self):
        var = float(self.line_Opump_eff.text())
        if var > 0.0 and var <= 1.0:
            main.default.Eff_po = var;
        else:
            self.line_Opump_eff.setText("0.3")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid oxidizer pump efficiency, try again.")
            msg.exec_()

    def checkFpumpEff(self):
        var = float(self.line_Fpump_eff.text())
        if var > 0.0 and var <= 1.0:
            main.default.Eff_pf = var;
        else:
            self.line_Fpump_eff.setText("0.3")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid fuel pump efficiency, try again.")
            msg.exec_()

    def checkTurbEff(self):
        var = float(self.line_turb_eff.text())
        if var > 0.0 and var <= 1.0:
            main.default.Eff_t = var;
        else:
            self.line_turb_eff.setText("0.5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid turbine efficiency, try again.")
            msg.exec_()

    def checkMechEff(self):
        var = float(self.line_mech_eff.text())
        if var > 0.0 and var <= 1.0:
            main.default.Eff_m = var;
        else:
            self.line_mech_eff.setText("0.95")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid mechanical efficiency, try again.")
            msg.exec_()

    def checkElePow(self):
        var = float(self.line_Ele_power.text())
        if var > 0.0 and var <= 1.0e10:
            main.default.Wmotor = var;
        else:
            self.line_Ele_power.setText("1.0e6")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid electric motor power, try again.")
            msg.exec_()

    def checkValveLoss(self):
        var = float(self.line_valve_losses.text())
        if var > 0.0 and var <= 1.0e6:
            main.default.v_loss = var;
        else:
            self.line_valve_losses.setText("1000")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid valve losses, try again.")
            msg.exec_()

    def checkDeMax(self):
        var = float(self.line_De_max.text())
        if var > 0.0 and var <= 100:
            main.default.De_max = var;
        else:
            self.line_De_max.setText("5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid maximum exit diameter, try again.")
            msg.exec_()

    def checkEpsMax(self):
        var = float(self.line_eps_max.text())
        if var > 0.0 and var <= 1000:
            main.default.Eps_max = var;
        else:
            self.line_eps_max.setText("300")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid maximum expansion ratio, try again.")
            msg.exec_()

    def checkConvAng(self):
        var = float(self.line_conv_ang.text())
        if var > 0.0 and var < 90:
            main.default.Theta_con = var;
        else:
            self.line_conv_ang.setText("60")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid convergent angle, try again.")
            msg.exec_()

    def checkDivAng(self):
        var = float(self.line_div_ang.text())
        if var > 0.0 and var < 90:
            main.default.Theta_conical = var;
        else:
            self.line_div_ang.setText("15")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid divergent angle, try again.")
            msg.exec_()

    def checkParAng0(self):
        var = float(self.line_par_0ang.text())
        if var > 0.0 and var < 90:
            main.default.Theta_bell = var;
        else:
            self.line_par_0ang.setText("55")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid parabola initial angle, try again.")
            msg.exec_()

    def checkParAngf(self):
        var = float(self.line_par_fang.text())
        if var > 0.0 and var < float(self.line_par_0ang.text()):
            main.default.TH_exit_bell = var;
        else:
            self.line_par_fang.setText("55")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid parabola final angle, try again.")
            msg.exec_()

    def checkCurvConv(self):
        var = float(self.line_curv_conv.text())
        if var > 0.0 and var < 4.0:
            main.default.R_u_ratio = var;
        else:
            self.line_curv_conv.setText("1.5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid throat curvature (convergent), try again.")
            msg.exec_()

    def checkCurvDiv(self):
        var = float(self.line_curv_div.text())
        if var > 0.0 and var < 2.0:
            main.default.R_u_bell = var;
        else:
            self.line_curv_div.setText("1.5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid throat curvature (divergent), try again.")
            msg.exec_()

    def checkSF(self):
        var = float(self.line_SF.text())
        if var > 1.0 and var < 20.0:
            main.default.SF = var;
        else:
            self.line_SF.setText("1.2")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid safety factor, try again.")
            msg.exec_()

    def checkCd(self):
        var = float(self.line_Cd.text())
        if var > 0.0 and var < 2.0:
            main.default.Cd = var;
        else:
            self.line_Cd.setText("0.7")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid discharge coefficient, try again.")
            msg.exec_()

    def checktIgn(self):
        var = float(self.line_tign.text())
        if var > 0.0 and var < 4000.0:
            main.default.ignburntime = var;
        else:
            self.line_tign.setText("4")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid ignitor burn time, try again.")
            msg.exec_()

# Main
if __name__ == '__main__':
    #generate app
    app = QApplication(sys.argv)
    mainwindow = MainWindow()
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(mainwindow)
    widget.setWindowTitle("LRE Design Tool - 0.1")

    #Show
    widget.show()
    widget.showMaximized()

    #close
    try:
        sys.exit(app.exec_())
    except:
        print("Exiting")
