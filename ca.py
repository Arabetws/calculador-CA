import math
w = 377

class sistema_trifasico:
    def __init__(self,vmax, L_mH, C_uF, R, angulo):
        self.L_mH = L_mH #Valor do Indutor
        self.C_uF = C_uF #Valor do Capacitor
        self.R = R       #Valor da resistencia
        self.angulo = angulo #angulo
        self.Vmax = vmax 

    def tensao_RMS_fonte(self):
        #tensão RMS da fonte 
        Vrms = self.Vmax / math.sqrt(2)
        return Vrms
    
    def mostrar_rms(self):
        #mostrar no terminal o valor do RMS
        Vrms = self.tensao_RMS_fonte()
        print(f"-"*10,"Tensão RMS da fonte", 10* "-")
        print(f"Tensão RMS: {Vrms:.2f} V")
    
    def tensao_fasorial_fonte(self):
        Vrms = self.tensao_RMS_fonte() #Pega o valoe RMS
        ang_rad = math.radians(self.angulo)

        #forma retangular
        real  = Vrms * math.cos(ang_rad) #parte real 
        imag  = Vrms *math.sin(ang_rad) #parte imaginaria
        fasor_retangular = complex(real, imag) #fasor retangular na forma complexa
        sinal = "+" if imag >= 0 else "-" #determina o sinal de acordo com o angulo

        #forma polar
        angulo_polar = self.angulo #pega o angulo digitafo
        sinal_polar = "+" if angulo_polar >= 0 else "-" 

        return fasor_retangular,sinal, angulo_polar
    
    def mostrar_fasorial_da_fonte(self):
        fasor, sinal, angulo = self.tensao_fasorial_fonte()
        print(f"-" * 10, "Tensão fasorial da fonte", 10 * "-")
        print(f"Fasor (forma retangular): {fasor.real:.2f}{sinal}j{abs(fasor.imag):.2f} [V]")
        print(f"Fasor (forma polar): {abs(fasor):.2f}.e{sinal}j{angulo} [V]")

    def reatancias_cap_ind(self):
        self.Xc = 1 / (w * self.C_uF * 1e-6)
        self.Xl = w * self.L_mH*1e-3

    def mostrar_reatancia(self):
        self.reatancias_cap_ind()
        print(f"-" * 10, "Reatancia capacitiva e indutiva", 10 * "-")
        print(f'Xc = {self.Xc:.2f}[Ω]')
        print(f'Xl = {self.Xl:.2f}[Ω]')
    
    def impedancia_equivalente(self):
        self.reatancias_cap_ind()
        Zeq = complex(self.R, self.Xl - self.Xc)
        real = Zeq.real
        imag = Zeq.imag
        sinal = "+" if imag >= 0 else "-"
        #forma polar
        modulo = abs(Zeq)
        angulo_rad = math.atan2(imag, real)
        angulo_graus = math.degrees(angulo_rad)
        sinal = "+" if angulo_graus>= 0 else "-"
        return Zeq, angulo_graus, real, sinal, imag, modulo, 
    
    def mostrar_impedancia_eq(self):
        Zeq,angulo, real, sinal, imag, modulo= self.impedancia_equivalente()
        print(f"-"*10, "Impedancia equivalente", 10 * "-")
        print(f'Forma Retangular = {real:.2f} {sinal} J {abs(imag):.2f}[Ω]')
        print(f'Forma Polar: {modulo:.2f}.e^{sinal}j{angulo:.2f}[Ω]')
    
    def corrente_fasorial(self):
        Vrms = self.tensao_RMS_fonte()
        _, _, angulo_fonte = self.tensao_fasorial_fonte()
        Zeq, angulo_impedancia, *_ = self.impedancia_equivalente()
            
        modulo_Z = abs(Zeq)
        corrente_modulo = Vrms / modulo_Z
        angulo_corrente = angulo_fonte - angulo_impedancia

        sinal = "+" if angulo_corrente >= 0 else "-"

        real = corrente_modulo *math.cos(angulo_corrente)
        img = corrente_modulo *math.sin(angulo_corrente)
        return corrente_modulo,angulo_corrente, sinal, real, img
    
    def mostrar_corrente_fasorial(self):
        corrente_modulo, angulo_corrente, sinal, real, img = self.corrente_fasorial()
        print(f"-"*10, "Corrente fasorial", 10* "-")
        print(f'i = {corrente_modulo:.2f}.e^{sinal}j{abs(angulo_corrente):.2f}[A]')
        print(f'i = {real:.2f}{sinal}j{img:.2f}[A]')

    def quedas_tensao(self):
        corrente_modulo, angulo_corrente, *_ = self.corrente_fasorial()
        self.impedancia_equivalente()
        #resistor
        Vr = self.R * corrente_modulo
        real_res = Vr *math.cos(angulo_corrente)
        img_res = Vr *math.sin(angulo_corrente)
        sinal = "+" if angulo_corrente >= 0 else "-"
        #capacitor
        Vc = self.Xc * corrente_modulo
        angulo_capacitor = - 90 +angulo_corrente
        ang_cap_rad = math.radians(angulo_capacitor)
        real_capacitivo =  Vc *math.cos(ang_cap_rad)
        img_cap =Vc *math.sin(ang_cap_rad)
        sinal_cap = "+" if img_cap >= 0 else "-"
        #indutor
        Vl = self.Xl * corrente_modulo
        angulo_indutor = 90 + angulo_corrente
        ang_ind_rad = math.radians(angulo_indutor)
        real_ind = Vl *math.cos(ang_ind_rad)
        img_ind = Vl *math.sin(ang_ind_rad)
        sinal_ind = "+" if img_ind >=0 else "-"
        return angulo_corrente, angulo_indutor, angulo_capacitor, real_res, real_capacitivo, real_ind, sinal, img_res, sinal_cap, img_ind, img_cap, sinal_ind, Vr, Vc, Vl

    def mostrar_queda_tensao(self):
        angulo_corrente, angulo_indutor, angulo_capacitor, real_res, real_capacitivo, real_ind, sinal, img_res, sinal_cap, img_ind, img_cap, sinal_ind, Vc, Vl, Vr = self.quedas_tensao() 
        print("-"*15, "Quedas de tensão", 15*"-")
        print("-" * 40)
        print(f"{'Componente':<15}{'Forma Polar':<20}{'Forma Retangular'}")
        print("-" * 40)
        print(f"{'Resistor':<15}{f'{Vr:.2f}.ej{angulo_corrente:.2f}':<20}{f'{real_res:.2f} {sinal}j{abs(img_res):.2f}'}")
        print(f"{'Capacitor':<15}{f'{Vc:.2f}.ej{angulo_capacitor:.2f}':<20}{f'{real_capacitivo:.2f} {sinal_cap}j{abs(img_cap):.2f}'}")
        print(f"{'Indutor':<15}{f'{Vl:.2f}.ej{angulo_indutor:.2f}':<20}{f'{real_ind:.2f} {sinal_ind}j{abs(img_ind):.2f}'}")
        print("-" * 40)
    def potencias_ativas(self):
        Vrms = self.tensao_RMS_fonte()
        corrente_modulo, angulo_corrente, *_ = self.corrente_fasorial()
        _, _, angulo_tensao = self.tensao_fasorial_fonte()

        # Diferença de fase em radianos
        delta_angulo_rad = math.radians(angulo_tensao - angulo_corrente)

        # Potências
        P = Vrms * corrente_modulo * math.cos(delta_angulo_rad)  # Ativa
        Q = Vrms * corrente_modulo * math.sin(delta_angulo_rad)  # Reativa
        S = Vrms * corrente_modulo                                # Aparente

        return Q,S,P

    def mostar_potencia_ativas(self):
        Q, P, S = self.potencias_ativas()
        print("-" * 10, "Potências", 10 * "-")
        print(f"P (ativa)     = {P:.2f} [W]")
        print(f"Q (reativa)   = {Q:.2f} [Var]")
        print(f"S (aparente)  = {S:.2f} [VA]")

vmax = float(input("Digite o valor de vmax(V): "))
L_mH = float(input("Digite o valor de L(mH): "))
C_uF = float(input("Digite o valor de C(µF): "))
R = float(input("Digite o valor de R(Ω): "))
angulo = float(input("Digite o valor do ângulo (graus): "))

# Cria o sistema trifásico com os valores fornecidos
sistema = sistema_trifasico(
    vmax=vmax,
    L_mH=L_mH,
    C_uF=C_uF,
    R=R,
    angulo=angulo,
)

# Executa os métodos
sistema.mostrar_rms()
sistema.mostrar_fasorial_da_fonte()
sistema.mostrar_reatancia()
sistema.mostrar_impedancia_eq()
sistema.mostrar_corrente_fasorial()
sistema.mostrar_queda_tensao()
sistema.mostar_potencia_ativas()

