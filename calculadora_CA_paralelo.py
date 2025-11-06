import math
w=377
class Corrente_AC_paralelo:
    def __init__(self,R1,R2, C_uF, L_mH, Vmax, angulo):
        self.R1 = R1
        self.R2 = R2
        self.C_uF = C_uF
        self.L_mH = L_mH
        self.Vmax = Vmax
        self.angulo = angulo #angulo de fase

    def tensao_rms(self):
        angulo_rad = math.radians(self.angulo)
        real = self.Vmax * math.cos(angulo_rad)
        imag = self.Vmax * math.sin(angulo_rad)
        sinal = "+" if imag >= 0 else "-"
        return real, imag, sinal

    def mostrar_tensao_rms(self):
        real, imag, sinal = self.tensao_rms()
        print(f'V={self.Vmax}.eJ{self.angulo}[v]')
        print(f'{real:.2f}{sinal}J{abs(imag):.2f}')
    def xl_e_xc(self):
        Xl = w * self.L_mH * 1e-3
        angulo_xl = 90
        Xc = 1 / (w * self.C_uF * 1e-6)
        angulo_xc = -90
        return Xl, Xc, angulo_xc, angulo_xl

    def impedancias(self):
        Xl, Xc, *_ = self.xl_e_xc()
        z1 = math.sqrt(self.R1**2 + Xl**2)
        z2 = math.sqrt(self.R2**2 + (-Xc)**2)
        angulo_z1_rad = math.atan(Xl / self.R1)
        angulo_z2_rad = math.atan(-Xc / self.R2)
        angulo_z1_graus = math.degrees(angulo_z1_rad)
        angulo_z2_graus = math.degrees(angulo_z2_rad)
        sinal_z1 = "+" if angulo_z1_graus >= 0 else "-"
        sinal_z2 = "+" if angulo_z2_graus >= 0 else "-"
        return z1, z2, angulo_z1_graus, angulo_z2_graus, sinal_z1, sinal_z2

    def mostrar_impedancia(self):
        z1, z2, angulo_z1, angulo_z2, sinal_Z1, sinal_z2 = self.impedancias()
        xl, xc, *_ = self.xl_e_xc()
        print(f'-' * 10, 'As impedâncias 𝑍1 de 𝑍̇2', 10 * '-')
        print(f'-' * 10, 'FORMA POLAR', 10 * '-')
        print(f'Z1 = {z1:.2f}.e{sinal_Z1}J{abs(angulo_z1):.2f}[Ω]')
        print(f'Z2 = {z2:.2f}.e{sinal_z2}J{abs(angulo_z2):.2f}[Ω]')
        print(f'-' * 10, 'FORMA RETANGULAR', 10 * '-')
        print(f'Z1 = {self.R1}+J{xl:.2f}[Ω]')
        print(f'Z2 = {self.R2}-J{xc:.2f}[Ω]')

    def correntes(self):
        z1, z2, angulo_z1, angulo_z2, *_ = self.impedancias()
        i1 = Vmax / z1
        i2 = Vmax / z2
        angulo_i1 = angulo_z1 + self. angulo
        angulo_i2 = angulo_z2 + self.angulo
        real_i1 = i1 *math.cos(angulo_i1)
        imag_i1 = i1 *math.sin(angulo_i1)

        real_i2 = i2 *math.cos(angulo_i2)
        imag_i2 = i2 *math.sin(angulo_i2)
        sinal_i1 = "+" if angulo_i1 >= 0 else "-"
        sinal_i2 = "+" if angulo_i2 >= 0 else "-"
        sinal_i1_ret = "+" if imag_i1 >=0 else "-"
        sinal_i2_ret = "+" if imag_i2 >=0 else "-"
        return i1, i2, angulo_i1, angulo_i2, real_i1, real_i2, imag_i1, imag_i2, sinal_i1, sinal_i2, sinal_i1_ret, sinal_i2_ret
    def mostrar_correntes(self):
        i1, i2, angulo_i1, angulo_i2, real_i1, real_i2, imag_i1, imag_i2, sinal_i1, sinal_i2, sinal_i1_ret, sinal_i2_ret = self.correntes()
        print(f'-' * 10, 'FORMA POLAR', 10 * '-')
        print(f'-' * 10, 'correntes 𝐼1 e 𝐼2', 10 * '-')
        print(f'i1 = {i1:.2f}.e{sinal_i1}J{angulo_i1:.2f}[A]')
        print(f'i2 = {i2:.2f}.e{sinal_i2}J{angulo_i2:.2f}[A]')
        print(f'-' * 10, 'FORMA RETANGULAR', 10 * '-')
        print(f'i1 = {real_i1:.2f}{sinal_i1_ret}J{abs(imag_i1):.2f}[A]')
        print(f'i2 = {real_i2:.2f}{sinal_i2_ret}J{abs(imag_i2):.2f}[A]')
    def quedas_de_tensao(self):
        i1, i2, angulo_i1, angulo_i2, *_ = self.correntes()
        xc, xl, *_ = self.xl_e_xc()
        # Tensões em forma polar
        Vr1 = self.R1 * i1
        Vr2 = self.R2 * i2
        Vc = xc * i2
        Vl = xl * i1

        # Ângulos (iguais aos das correntes)
        vr1_angulo = angulo_i1
        vr2_angulo = angulo_i2
        vc_angulo = angulo_i2
        vl_angulo = angulo_i1

        # Sinais para exibição
        sinal_R1_pol = "+" if vr1_angulo >= 0 else "-"
        sinal_R2_pol = "+" if vr2_angulo >= 0 else "-"
        sinal_vc_pol = "+" if vc_angulo >= 0 else "-"
        sinal_vl_pol = "+" if vl_angulo >= 0 else "-"

        # Conversão para radianos
        vr1_rad = math.radians(vr1_angulo)
        vr2_rad = math.radians(vr2_angulo)
        vc_rad = math.radians(vc_angulo)
        vl_rad = math.radians(vl_angulo)

        # Forma retangular
        real_r1 = Vr1 * math.cos(vr1_rad)
        imag_r1 = Vr1 * math.sin(vr1_rad)
        sinal_r1_ret = "+" if imag_r1 >= 0 else "-"

        real_r2 = Vr2 * math.cos(vr2_rad)
        imag_r2 = Vr2 * math.sin(vr2_rad)
        sinal_r2_ret = "+" if imag_r2 >= 0 else "-"

        real_vc = Vc * math.cos(vc_rad)
        imag_vc = Vc * math.sin(vc_rad)
        sinal_vc_ret = "+" if imag_vc >= 0 else "-"

        real_vl = Vl * math.cos(vl_rad)
        imag_vl = Vl * math.sin(vl_rad)
        sinal_vl_ret = "+" if imag_vl >= 0 else "-"

        return {
            "Vr1_polar": (Vr1, vr1_angulo, sinal_R1_pol),
            "Vr2_polar": (Vr2, vr2_angulo, sinal_R2_pol),
            "Vc_polar": (Vc, vc_angulo, sinal_vc_pol),
            "Vl_polar": (Vl, vl_angulo, sinal_vl_pol),

            "Vr1_retangular": (real_r1, imag_r1, sinal_r1_ret),
            "Vr2_retangular": (real_r2, imag_r2, sinal_r2_ret),
            "Vc_retangular": (real_vc, imag_vc, sinal_vc_ret),
            "Vl_retangular": (real_vl, imag_vl, sinal_vl_ret)
        }
    def mostrar_queda_tensao(self):
        tensoes = self.quedas_de_tensao()

        print("-" * 40)
        print("TABELA DE QUEDAS DE TENSÃO")
        print("-" * 40)
        print(f"{'Elemento':<12} {'Polar':<20} {'Retangular':<20}")
        print("-" * 40)

        # Resistor R1
        vr1, angulo_r1, sinal_r1_pol = tensoes["Vr1_polar"]
        real_r1, imag_r1, sinal_r1_ret = tensoes["Vr1_retangular"]
        print(f"{'R1':<12} {vr1:.2f}.e{sinal_r1_pol}J{abs(angulo_r1):.2f} {real_r1:.2f}{sinal_r1_ret}J{abs(imag_r1):.2f}")

        # Resistor R2
        vr2, angulo_r2, sinal_r2_pol = tensoes["Vr2_polar"]
        real_r2, imag_r2, sinal_r2_ret = tensoes["Vr2_retangular"]
        print(f"{'R2':<12} {vr2:.2f}.e{sinal_r2_pol}J{abs(angulo_r2):.2f}  {real_r2:.2f}{sinal_r2_ret}J{abs(imag_r2):.2f}")

        # Capacitor
        vc, angulo_c, sinal_c_pol = tensoes["Vc_polar"]
        real_c, imag_c, sinal_c_ret = tensoes["Vc_retangular"]
        print(f"{'Capacitor':<12} {vc:.2f}.e{sinal_c_pol}J{abs(angulo_c):.2f}  {real_c:.2f}{sinal_c_ret}J{abs(imag_c):.2f}")

        # Indutor
        vl, angulo_l, sinal_l_pol = tensoes["Vl_polar"]
        real_l, imag_l, sinal_l_ret = tensoes["Vl_retangular"]
        print(f"{'Indutor':<12} {vl:.2f}.e{sinal_l_pol}J{abs(angulo_l):.2f}  {real_l:.2f}{sinal_l_ret}J{abs(imag_l):.2f}")

        print("-" * 40)



Vmax = float(input("Digite o valor de vmax(V): "))
L_mH = float(input("Digite o valor de L(mH): "))
C_uF = float(input("Digite o valor de C(µF): "))
R1= float(input("Digite o valor do primeiro resistor(Ω): "))
R2 =float(input("Digite o valor do segundo resistor(Ω): "))
angulo = float(input("Digite o valor do ângulo (graus): "))

sistema = Corrente_AC_paralelo(
    R1 = R1,
    R2 = R2,
    C_uF = C_uF,
    L_mH = L_mH,
    Vmax = Vmax,
    angulo = angulo
)

sistema.mostrar_tensao_rms()
sistema.mostrar_impedancia()
sistema.mostrar_correntes()
sistema.mostrar_queda_tensao()
