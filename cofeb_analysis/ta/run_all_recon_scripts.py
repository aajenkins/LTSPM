# @Author: Jenkins Alec <alec>
# @Date:   2017-08-03T19:15:33-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-10T10:16:47-07:00



import vector_m_reconstruction
import calc_DW
import stray_field
import plot_stray_field

def run_all(scannum, fitHelicity=False, h=45):
    if (fitHelicity):
        helicity_list=[0,90,180,h]
    else:
        helicity_list=[0,90,180]

    print('calc stray vector field and Mz')
    vector_m_reconstruction.vector_m_reconstruction(scannum)
    print('calc magnetization')
    calc_DW.calc_DW(scannum, helicities=helicity_list)
    stray_field.stray_field(scannum, helicities=helicity_list)
    plot_stray_field.plot_stray_field(scannum, helicities=helicity_list)


if __name__ == "__main__":
    import sys
    if (len(sys.argv) == 2):
        run_all(int(sys.argv[1]))
    elif (len(sys.argv) == 4):
        run_all(int(sys.argv[1]), fitHelicity=eval(sys.argv[2]), h=eval(sys.argv[3]))
    else:
        print('enter scan number or scan number and both fitHelicity bool and helicity number')
