import math
import matplotlib.pyplot as plt
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
        corrente_modulo, angulo_corrente_rad, *_ = self.corrente_fasorial()
        self.impedancia_equivalente()
        
        # --- Resistor (VR) ---
        Vr_mag = self.R * corrente_modulo
        # VR está em fase com a corrente (ângulo = angulo_corrente_rad)
        real_res = Vr_mag * math.cos(angulo_corrente_rad)
        img_res = Vr_mag * math.sin(angulo_corrente_rad)
        
        # --- Capacitor (VC) ---
        Vc_mag = self.Xc * corrente_modulo
        # VC atrasa a corrente em 90 graus (ângulo = angulo_corrente_rad - 90 graus)
        angulo_capacitor_rad = angulo_corrente_rad - math.radians(90)
        real_capacitivo = Vc_mag * math.cos(angulo_capacitor_rad)
        img_cap = Vc_mag * math.sin(angulo_capacitor_rad)
        
        # --- Indutor (VL) ---
        Vl_mag = self.Xl * corrente_modulo
        # VL adianta a corrente em 90 graus (ângulo = angulo_corrente_rad + 90 graus)
        angulo_indutor_rad = angulo_corrente_rad + math.radians(90)
        real_ind = Vl_mag * math.cos(angulo_indutor_rad)
        img_ind = Vl_mag * math.sin(angulo_indutor_rad)
        
        # Retornando os componentes reais e imaginários, magnitudes e ângulos em graus
        return real_res, img_res, real_capacitivo, img_cap, real_ind, img_ind, Vr_mag, Vc_mag, Vl_mag, angulo_corrente_rad, angulo_indutor_rad, angulo_capacitor_rad
    
    def diagrama(self):
        # Desempacota os resultados de quedas_tensao
        real_res, img_res, real_capacitivo, img_cap, real_ind, img_ind, Vr_mag, Vc_mag, Vl_mag, I_ang_deg, Vl_ang_deg, Vc_ang_deg = self.quedas_tensao()

        # 1. Configuração da figura com 3 subplots
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle('Diagramas Fasoriais das Quedas de Tensão por Componente', fontsize=16)

        # Determinar o limite máximo para os eixos (para manter a proporção)
        max_val = max(Vr_mag, Vc_mag, Vl_mag) * 1.2 # 20% de margem
        lim = max(10, round(max_val / 10) * 10) # Arredonda para o 10 mais próximo, mínimo de 10

        # --- Função auxiliar para desenhar o diagrama ---
        def desenhar_diagrama(ax, real, imag, mag, ang_deg, titulo, cor):
            # Configurações básicas
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_xlabel('Parte Real (V)')
            ax.set_ylabel('Parte Imaginária (V)')
            ax.set_title(titulo)
            ax.grid(True)
            ax.axhline(0, color='black', linewidth=0.5)
            ax.axvline(0, color='black', linewidth=0.5)
            ax.set_aspect('equal', adjustable='box')

            # Desenha o vetor (fasor)
            ax.quiver(0, 0, real, imag, angles='xy', scale_units='xy', scale=1, color=cor, label=titulo, zorder=5)

            # Anotação do ângulo (em graus)
            ang_text = f'{ang_deg:.2f}°'
            
            # Desenha o arco do ângulo (opcional, mas ajuda a visualizar)
            start_ang = 0
            end_ang = ang_deg
            if ang_deg < 0:
                start_ang = ang_deg
                end_ang = 0
            
            # Desenha o arco
            if abs(ang_deg) > 0.1:
                arc_radius = lim * 0.15
                arc = plt.matplotlib.patches.Arc((0, 0), 2*arc_radius, 2*arc_radius, 
                                                angle=0, theta1=start_ang, theta2=end_ang, 
                                                color='gray', linestyle=':', linewidth=1)
                ax.add_patch(arc)
                
                # Posição do texto do ângulo no meio do arco
                mid_ang_rad = math.radians((start_ang + end_ang) / 2)
                text_x = arc_radius * 1.1 * math.cos(mid_ang_rad)
                text_y = arc_radius * 1.1 * math.sin(mid_ang_rad)
                
                ax.text(text_x, text_y, ang_text, 
                        fontsize=9, color='black', ha='center', va='center', 
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3'))


            # Anotação da magnitude (próximo ao vetor)
            ax.text(real * 0.8, imag * 0.8, f'{mag:.2f} V', 
                    fontsize=10, color=cor, ha='center', va='center', 
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3'))
            
            # Projeções nos eixos (linhas tracejadas)
            ax.plot([real, real], [0, imag], 'k--', linewidth=0.5) # Projeção Imaginária
            ax.plot([0, real], [imag, imag], 'k--', linewidth=0.5) # Projeção Real
            
            # Anotações das componentes (ajustadas para não sobrepor os eixos)
            ax.text(real, -lim * 0.05, f'{real:.2f}', fontsize=8, ha='center', va='top')
            ax.text(-lim * 0.05, imag, f'{imag:.2f}', fontsize=8, ha='right', va='center')


        # --- 2. Diagrama do Resistor (VR) ---
        desenhar_diagrama(axes[0], real_res, img_res, Vr_mag, I_ang_deg, 'Resistor (VR)', 'red')

        # --- 3. Diagrama do Capacitor (VC) ---
        desenhar_diagrama(axes[1], real_capacitivo, img_cap, Vc_mag, Vc_ang_deg, 'Capacitor (VC)', 'blue')

        # --- 4. Diagrama do Indutor (VL) ---
        desenhar_diagrama(axes[2], real_ind, img_ind, Vl_mag, Vl_ang_deg, 'Indutor (VL)', 'green')

        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Ajusta o layout para evitar sobreposição de títulos
        plt.show()

    def mostrar_queda_tensao(self):
        # Desempacota os novos retornos de quedas_tensao
        # real_res, img_res, real_capacitivo, img_cap, real_ind, img_ind, 
        # Vr, Vc, Vl, 
        # angulo_corrente, angulo_indutor, angulo_capacitor
        real_res, img_res, real_capacitivo, img_cap, real_ind, img_ind, \
        Vr, Vc, Vl, \
        angulo_corrente, angulo_indutor, angulo_capacitor = self.quedas_tensao()

        # Recalcula os sinais para a forma retangular
        sinal = "+" if img_res >= 0 else "-"
        sinal_cap = "+" if img_cap >= 0 else "-"
        sinal_ind = "+" if img_ind >= 0 else "-"

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
sistema.diagrama()
