# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 10:14:50 2022

@author: psypple
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xlwt

'''定义初始参数'''
table=pd.read_excel("structure.xlsx")
λ=table['波长（nm）'].values.tolist()[0]*1E-9
k0=2*np.pi/λ+0j
b=1E-6 #迭代用参数，决定初始步长，太大会导致跳出波导解
δ=1E-10 #迭代用参数
ε=1E-16 #迭代用参数
iter_max=1000 #设置一个最大迭代步数
iter_max2=100 #设置二级迭代，即漏模不收敛时，增加有源区增益的最大迭代次数
iter_max2_reduce=500 #设置增益倍增后缩减增益的最大迭代次数
net_gain=1E-6 #计算结果净增益值上限
δ_gain_init=100*100*λ/(4*np.pi) #收敛后，迭代增益使得净增益为0的初始步长
active_index_list=table['是否为有源区（是为1）'].values.tolist() #记录是否为有源区
active_index=[]
for i in range(len(active_index_list)):
    if active_index_list[i]==1:
        active_index.append(i)

'''输入折射率和厚度'''
n_real=table['折射率实部'].values.tolist()
n_imag=table['吸收（/cm）'].values.tolist()
for i in range(len(n_imag)):
    n_imag[i]=n_imag[i]*100*λ/(4*np.pi)
global n,d_m
n=[1+0j] #折射率，注意这个数组index=0时为空气，index最大时为衬底
for i in range(len(n_imag)):
    n.append(n_real[i]+n_imag[i]*1j)
d=table['层厚（nm）'].values.tolist() #每层厚度（从上往下）
del d[-1] #删除多余的nan项，这一项是因为无限衬底的厚度没有输入值
d_m=d #约化为国际标准单位m
for i in range (len(d_m)):
    d_m[i]=d_m[i]/1E9

'''三层波导试探解'''
d_waveguide=table['波导厚度（nm）'].values.tolist()[0]*1E-9 #波导总厚度
n1=table['cladding折射率实部'].values.tolist()[0] #cladding折射率，注意只要实部
n2=table['波导折射率实部'].values.tolist()[0] #波导折射率，注意只要实部
R_square=k0.real**2*d_waveguide**2/4*(n2**2-n1**2)
x_0=[]
δy_even=[]
δy_odd=[]
for i in range(10000):
    x_0.append(i*np.sqrt(R_square)/10000)
    if abs(np.sin(x_0[i]))<1E-10: #防止tan导致的发散
        δy_even.append(abs(np.sqrt(R_square-x_0[i]**2)-x_0[i]*1E-10))
        δy_odd.append(abs(np.sqrt(R_square-x_0[i]**2)+x_0[i]/1E-10))
    elif abs(np.cos(x_0[i]))<1E-10: #防止tan导致的发散
        δy_even.append(abs(np.sqrt(R_square-x_0[i]**2)-x_0[i]*1E10))
        δy_odd.append(abs(np.sqrt(R_square-x_0[i]**2)+x_0[i]/1E10))
    else:
        δy_even.append(abs(np.sqrt(R_square-x_0[i]**2)-x_0[i]*np.tan(x_0[i])))
        δy_odd.append(abs(np.sqrt(R_square-x_0[i]**2)+x_0[i]/np.tan(x_0[i])))
n0=[]
for i in range(9998):
    if δy_even[i+1]<δy_even[i] and δy_even[i+1]<δy_even[i+2]:
        n0.append(np.sqrt(n2**2*k0**2-4*x_0[i+1]**2/d_waveguide**2))
    if δy_odd[i+1]<δy_odd[i] and δy_odd[i+1]<δy_odd[i+2]:
        n0.append(np.sqrt(n2**2*k0**2-4*x_0[i+1]**2/d_waveguide**2))
for i in range(len(n0)):
    n0[i]=n0[i]/k0 #这里n0记录了所有的有效折射率试探解的解析解
print('Trial solution(s) for neff:')
for i in range(len(n0)):
    print('%f'%n0[i].real)
print ('substrate refractive index: %f' %n_real[-1])
print()

'''传递矩阵'''
def mode_profile(β):
    layer_num=len(d)
    γ=[] #γ即生长方向的各层中光场的复波矢
    for i in range(layer_num+2):
        γ.append(np.sqrt(β**2-k0**2*n[i]**2))
    AB=[np.array([1,0])] #A和B为光场通解的系数，AB列表记录了每一层（从0层空气开始）的每一组系数[A,B]
    for i in range(layer_num): #迭代layer_num次后，AB的最后一个组[Az-1,Bz-1]其实是倒数第二层，还没有得到无限衬底的系数
        Ti=np.array([[(1+γ[i]/γ[i+1])*np.exp(γ[i+1]*d_m[i])/2 , (1-γ[i]/γ[i+1])*np.exp(γ[i+1]*d_m[i])/2],
                     [(1-γ[i]/γ[i+1])*np.exp(-γ[i+1]*d_m[i])/2 , (1+γ[i]/γ[i+1])*np.exp(-γ[i+1]*d_m[i])/2]])
        AB_i1=np.dot(Ti,AB[i])
        AB.append(AB_i1)
    Tzmin1=np.array([[(1+γ[layer_num]/γ[layer_num+1])/2 , (1-γ[layer_num]/γ[layer_num+1])/2],
                     [(1-γ[layer_num]/γ[layer_num+1])/2 , (1+γ[layer_num]/γ[layer_num+1])/2]])
    AB.append(np.dot(Tzmin1,AB[-1])) #这里加入的是无限衬底层的系数，但是参考点与前面不同，计算光场时要单独处理
    return AB[layer_num+1][0],AB,γ #输出的第一项是Az，由于设置了A0=0，因此Az=t11用于迭代求解；后两项用于计算光场

'''迭代'''
def downhill(β_init):
    converge_state=0 #记录迭代是否收敛的值，不收敛会增加有源区增益
    β_iter=[β_init] #记录迭代中的β，最后一项是最终符合收敛条件的β，即解
    t11=abs(mode_profile(β_init)[0])
    t11_0=[t11] #记录迭代中的t11，最后一项是最终符合收敛条件的t11<δ
    δβ=b*t11
    crimdef=0 #这个参数叫做前科值，即上次的δβ过大造成t11迭代没有减小留下记录，下次迭代若成果缩小则消除记录，但当此不减半δβ
    count=0
    while t11>=δ and count<iter_max:
        count=count+1
        β_temp=[] #方便对照四个t11'中哪个是下降最快的解
        t11_1=[] #记录四个试探解的t11的模
        t11_11=mode_profile(β_iter[-1]+δβ)[0]
        t11_12=mode_profile(β_iter[-1]-δβ)[0]
        t11_13=mode_profile(β_iter[-1]+δβ*1j)[0]
        t11_14=mode_profile(β_iter[-1]-δβ*1j)[0]
        t11_1.append(abs(t11_11))
        t11_1.append(abs(t11_12))
        t11_1.append(abs(t11_13))
        t11_1.append(abs(t11_14))
        β_temp.append(β_iter[-1]+δβ)
        β_temp.append(β_iter[-1]-δβ)
        β_temp.append(β_iter[-1]+1j*δβ)
        β_temp.append(β_iter[-1]-1j*δβ)
        if min(t11_1)<t11:
            β_iter.append(β_temp[t11_1.index(min(t11_1))])
            t11=t11_1[t11_1.index(min(t11_1))]
            t11_0.append(t11)
            if crimdef==0:
                δβ=1.1*δβ
            crimdef=0
        else:
            crimdef=1
            if not δβ<ε:
                δβ=0.5*δβ
    if t11>=δ:
        converge_state=0
    else:
        converge_state=1
    for i in range (len(β_iter)):
        β_iter[i]=β_iter[i]/k0
    return β_iter,t11_0,converge_state #此时输出的是有效折射率

'''主控制脚本'''
def main(neff_0,q,plot): #主程序要求输入迭代初始的有效折射率，注意初始值非常关键！特别是存在高阶模时，输入的第二项q是为了记录解的数量保存图片用
    q=q+1
    matrix=downhill(neff_0*k0) #导出迭代步骤数据
    β_final=matrix[0] #迭代过程中有效折射率中间数值
    t11=matrix[1] #迭代过程中t11中间数值
    converge_state=matrix[2] #传递downhill函数中的收敛性状态值
    main_ABγ=mode_profile(β_final[-1]*k0)
    iter_num=[]
    for i in range(len(β_final)):
        iter_num.append(i)
    # plt.figure()
    # plt.plot(iter_num,β_final)
    # plt.xlabel('iter number')
    # plt.ylabel('neff_iter')
    # plt.show()
    if plot==1:
        plt.figure()
        plt.title('Solution %d' %q)
        plt.plot(iter_num,t11)
        plt.xlabel('iter number')
        plt.ylabel('t11_iter')
        plt.yscale('log')
        # plt.savefig('Convergence for Solution %d.png' %q)
        plt.show()
    return main_ABγ,β_final[-1],converge_state

'''主求解程序'''
def soluter():
    δ_gain=δ_gain_init
    for q in range(len(n0)):
        q1=q+1
        ABγ=main(n0[q],q,0) #把解系数导出
        AB=ABγ[0][1]
        neff_imag=ABγ[1].imag
        converge_state=ABγ[2] #传递收敛状态值
        count=0
        noleak=converge_state #此值为1时跳过增益迭代，因为没有漏模
        while converge_state==0 and count<iter_max2: #不收敛时开始增加增益
            count=count+1
            for i in active_index:
                n[i+1]=n[i+1]-δ_gain*1j #首先使得增益非0，然后进行倍增直到收敛
            ABγ=main(n0[q],q,0) #重新求解
            AB=ABγ[0][1]
            neff_imag=ABγ[1].imag
            converge_state=ABγ[2]
            δ_gain=1.5*δ_gain
        if converge_state==0:
            print ('Failed to convert.')
        else:
            print ('Converting succeeds.')
            ABγ_temp=ABγ #临时记录收敛成功的解，如果增益迭代失败则引用这个解
            AB_temp=AB
            neff_imag_temp=neff_imag
            if noleak==0:
                δ_gain=δ_gain_init
                count=0
                crimedef=0
                while count<iter_max2_reduce and abs(neff_imag)>net_gain: #downhill迭代增益使得净增益为0
                    compare_neffimag=[] #分别对有源区折射率虚部加和减，然后比较哪样使得净增益减小
                    compare_converge=[]
                    for i in active_index: #先增加
                        n[i+1]=n[i+1]+δ_gain*1j
                    ABγ_plus=main(n0[q],q,0) #重新求解
                    AB_plus=ABγ_plus[0][1]
                    neff_imag_plus=ABγ_plus[1].imag
                    converge_state_plus=ABγ_plus[2]
                    compare_neffimag.append(abs(neff_imag_plus))
                    compare_converge.append(converge_state_plus)
                    for i in active_index: #先增加
                        n[i+1]=n[i+1]-2*δ_gain*1j
                    ABγ_minus=main(n0[q],q,0) #重新求解
                    AB_minus=ABγ_minus[0][1]
                    neff_imag_minus=ABγ_minus[1].imag
                    converge_state_minus=ABγ_minus[2]
                    compare_neffimag.append(abs(neff_imag_minus))
                    compare_converge.append(converge_state_minus)
                    if min(compare_converge)==1 and min(compare_neffimag)<abs(neff_imag):
                        if compare_neffimag[0]<compare_neffimag[1]:
                            ABγ=ABγ_plus
                            AB=AB_plus
                            neff_imag=neff_imag_plus
                            converge_state=converge_state_plus
                            for i in active_index: #先增加
                                n[i+1]=n[i+1]+2*δ_gain*1j
                            if crimedef==0:
                                δ_gain=1.2*δ_gain
                            else:
                                crimedef=0
                        else:
                            ABγ=ABγ_minus
                            AB=AB_minus
                            neff_imag=neff_imag_minus
                            converge_state=converge_state_minus
                            if crimedef==0:
                                δ_gain=1.2*δ_gain
                            else:
                                crimedef=0
                    else:
                        crimedef=1
                        for i in active_index:
                            n[i+1]=n[i+1]+δ_gain*1j
                        δ_gain=0.5*δ_gain
                    count=count+1
                if converge_state==0:
                    print ('Failed to achieve a steady state.')
                    ABγ=ABγ_temp #重新引用临时解（收敛）
                    AB=AB_temp
                    neff_imag=neff_imag_temp
        A=[]
        B=[]
        γ=ABγ[0][2]
        for i in range(len(n)):
            A.append(AB[i][0])
            B.append(AB[i][1])
        solut_num=q+1
        print ('Solution %d:'%solut_num)
        # print ('Az=%e + %e i' %(A[-1].real,A[-1].imag)) #输出t11的值判断精确度
        # print ('Bz=%e + %e i' %(B[-1].real,B[-1].imag)) #无限衬底开始处电场强度
        print ('neff=%f' %ABγ[1].real)
        
        '''绘制光场'''
        t=[0] #设置边界点，计算光场
        for i in range(len(d)):
            t.append(t[i]+d[i])
        x=[] #x是坐标，用于绘制光场，单位nm
        E=[] #E是光场
        n_plot=[] #画光场轮廓
        integral_active=[] #计算限制因子分子的积分
        for i in range(len(d)):
            integral_active.append(0)
        integral=0 #计算限制因子分母的积分
        for i in range(200): #上空气层光场
            x.append(1E-9*(-199+i))
            E.append(np.exp(γ[0]*(i-199)*1E-9))
            n_plot.append(n[0].real)
        for i in range(len(d)): #计算有限厚度层光场
            for k in range(500):
                x.append(t[i]+k*d[i]/500)
                Elec_field=A[i+1]*np.exp(γ[i+1]*(k*d[i]/500-d[i]))+B[i+1]*np.exp(-γ[i+1]*(k*d[i]/500-d[i]))
                E.append(Elec_field)
                n_plot.append(n[i+1].real)
                if len(active_index)>=1: #设置存在有源区时计算限制因子
                    if i in active_index:
                        integral_active[i]=integral_active[i]+d[i]/500*abs(Elec_field)**2
        if γ[-1].real>1E-10:
            x_submax=10*np.log(10)/abs(γ[-1].real) #设置衬底计算最远距离防止光场发散
        else:
            x_submax=3.1E-6
        if x_submax>3E-6:
            x_submax=3E-6
        E_structure=E #记录没有考虑无限衬底时的光场，因为计算限制因子不能加衬底
        I_modal_structure=[] #计算限制因子用，不加入衬底的光场强度
        for i in range(len(E_structure)):
            I_modal_structure.append(abs(E_structure[i])**2)
        for i in range(len(E_structure)-1):
            integral=integral+I_modal_structure[i]*(x[i+1]-x[i])
        for i in range(1000): #计算无限衬底光场
            x.append(t[-1]+i*x_submax/1000)
            E.append(A[-1]*np.exp(γ[-1]*i*x_submax/1000)+B[-1]*np.exp(-γ[-1]*i*x_submax/1000))
            n_plot.append(n[-1].real)
        E_real=[]
        E_imag=[]
        I_modal=[] #复光场的模
        E_phase=[] #复光场的相位
        for i in range(len(E)):
            E_real.append(E[i].real)
            E_imag.append(E[i].imag)
            I_modal.append(abs(E[i])**2)
            E_phase.append(E[i].imag/abs(E[i])+1.5) #注意+1.5是为了画图好看
        for i in range(len(x)):
            x[i]=x[i]*1E6
        Γ0=[]
        for i in range(len(d)):
            if i in active_index:
                Γ0.append(integral_active[i]/integral)
            else:
                Γ0.append(0)
        Γ=sum(Γ0)
        α=0
        for i in active_index:
            α=α-n[i+1].imag*(4*np.pi)/(100*λ)*Γ0[i]
        print ('Γ=%f' %Γ)
        print ('Total modal gain=%f /cm' %α)
        print ()
        fig, ax1=plt.subplots(1,1)
        ax2=ax1.twinx()
        plt.title('Solution %d, Γ=%f' %(q1,Γ))
        ax1.plot(x,I_modal,color='red')
        ax1.set_xlabel('Distance (μm)')
        ax1.set_ylabel('Intensity (arb.unit)',color='red')
        ax1.tick_params(axis='y',colors='red')
        ax2.plot(x,n_plot,color='blue')
        if ABγ[1].real<n_real[-1]: #仅仅存在漏模时画相位图
            ax2.plot(x,E_phase,color='red',linestyle='dashed')
        ax2.set_ylabel('refractive index',color='blue')
        ax2.set_ylim(0,max(n_plot)*1.5)
        ax2.tick_params(axis='y',colors='blue')
        ax2.spines['left'].set_color('red')
        ax2.spines['right'].set_color('blue')
        plt.savefig('Modal Intensity for Solution %d.png' %q1)
        plt.show()
        
        '''打印结果'''
        f = xlwt.Workbook() #创建工作薄
        sheet1 = f.add_sheet(u'sheet1',cell_overwrite_ok=True) #创建sheet
        list = x
        list.insert(0,'distance (μm)')
        list.insert(0,'neff')
        list.insert(0,'Γ')
        list.insert(0,'total modal gain (/cm)')
        j = 0
        for i in list:
            sheet1.write(j,0,i) #循环写入
            j=j+1
        list = I_modal
        list.insert(0,'Intensity (arb.unit)')
        list.insert(0,'%f' %ABγ[1].real)
        list.insert(0,'%f' %Γ)
        list.insert(0,'%f' %α)
        j = 0
        for i in list:
            sheet1.write(j,1,i) #循环写入
            j=j+1
        list = E_real
        list.insert(0,'Re(E) (arb.unit)')
        list.insert(0,'')
        list.insert(0,'')
        list.insert(0,'')
        j = 0
        for i in list:
            sheet1.write(j,2,i) #循环写入
            j=j+1
        list = E_imag
        list.insert(0,'Im(E) (arb.unit)')
        list.insert(0,'')
        list.insert(0,'')
        list.insert(0,'')
        j = 0
        for i in list:
            sheet1.write(j,3,i) #循环写入
            j=j+1
        f.save('solution%d.xls' %q1)#保存文件
soluter()
input('Press any button to exit.')