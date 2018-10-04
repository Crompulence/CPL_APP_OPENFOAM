import numpy as np

class StokesAnalytical:

	def __init__(self, rhof=1.00):

		self.rhof = rhof
		self.read_lammps_inputfiles()

	def read_lammps_inputfiles(self):
		"""

		Extract the MD timestep from the input file (single.in). Extract the
		particle diameter, density, initial position (in the y-direction) and
		initial velocity from the data file (single.lj). Note the file name
		and their location in the 'lammps' sub-directory is hard-wired.

		"""
		with open('lammps/single.in', 'r') as fObj:
			for l in fObj:
				if 'timestep' in l:
					dt = float(l.split()[1])
				if 'gravity' in l:
					g = -float(l.split()[4])
				if 'viscosity equal' in l:
					mu = float(l.split()[3])

		with open('lammps/single.lj', 'r') as fObj:
			lines = filter(None, (line.rstrip() for line in fObj))   
			for i, l in enumerate(lines):
				if l == 'Atoms':
					data = lines[i+1].split()
					dp = float(data[2])
					rhop = float(data[3])
					x0 = float(data[5])
					
				if l == 'Velocities':
					data = lines[i+1].split()
					v0 = float(data[5])

		self.x0 = x0
		self.v0 = v0
		self.dp = dp
		self.rhop = rhop
		self.dt = dt
		self.g = g
		self.mu = mu

	def read_thermo_output(self):
		"""
		
		Extract the time, velocity and displacement of the particle during the
		simulation from the 'thermo_output.txt' file in the 'lammps'
		directory. Note that this will not include the initial position and
		velocity, which must be inserted into the array.

		"""

		data = np.loadtxt('lammps/thermo_output.txt', skiprows=1)
		t = data[:,1]
		xp = data[:,3]
		vp = data[:,6]
		t = np.insert(t, 0, 0)
		xp = np.insert(xp, 0, self.x0)
		vp = np.insert(vp, 0, self.v0)

		self.t = t 
		self.xp = xp
		self.vp = vp

	def get_analytical_solution(self):
		"""

		Analytical displacement and velocity profile for single particle
		falling under gravity in a stationary fluid, along with the terminal
		velocity.

		"""

		self.read_thermo_output()

		mp = self.rhop*(np.pi/6)*self.dp**3
		mf = self.rhof*(np.pi/6)*self.dp**3
		kd = 3*np.pi*self.mu*self.dp

		vel = ((mp - mf)*self.g/kd)*(1 - np.exp(-kd*self.t/mp)) + self.v0*np.exp(-kd*self.t/mp)
		disp = ((mp - mf)*self.g/kd)*(self.t - (mp/kd)*(1 - np.exp(-kd*self.t/mp))) + (self.v0*mp/kd)*(1 - np.exp(-kd*self.t/mp)) + self.x0

		vel_term = (mp - mf)*self.g/kd

		self.disp = disp
		self.vel = vel 
		self.vel_term = vel_term

	def test_displacement(self, tol_rel=0.01, tol_abs=1.e-6, cutoff=1.e-3):
		"""

		Test displacement profile with analytical solution. If the analytical
		solution has a value less than the cutoff value, then perform a
		absolute error analysis, otherwise perform relative error analysis.

		"""

		xp = self.xp
		disp = self.disp		
		for idx in range(len(disp)):
			if abs(disp[idx]) < cutoff:
				assert(abs(disp[idx] - xp[idx]) <= tol_abs)
			else:
				assert(abs((disp[idx] - xp[idx])/disp[idx]) <= tol_rel)

	def test_velocity(self, tol_rel=0.01, tol_abs=1.e-6, cutoff=1.e-3):
		"""

		Test displacement profile with analytical solution. If the analytical
		solution has a value less than the cutoff value, then perform a
		absolute error analysis, otherwise perform relative error analysis.

		"""

		vp = self.vp
		vel = self.vel		
		for idx in range(len(vel)):
			if abs(vel[idx]) < cutoff:
				assert(abs(vel[idx] - vp[idx]) <= tol_abs)
			else:
				assert(abs((vel[idx] - vp[idx])/vel[idx]) <= tol_rel)

	def test_terminal(self, tol_rel=0.01):
		"""

		Test that terminal velocity has been obtained with an relative error
		tolerance.

		"""

		vp = self.vp[-1]
		err_rel = abs((vp - self.vel_term)/self.vel_term) <= tol_rel
		assert(err_rel)


	def plot_displacement_velocity(self):
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt

		plt.plot(self.t, self.vp, 'r-', linewidth=3.0)
		plt.plot(self.t, self.vel, 'k-')
		plt.plot(self.t, self.vel_term*np.ones_like(self.t), 'k--')
		plt.xlabel('Time (s)')
		plt.ylabel('Velocity (cm/s)')
		plt.legend(('Numerical', 'Analytical', 'Terminal'))
		plt.savefig('fig_velocity.png')
		plt.close()

		plt.plot(self.t, self.xp, 'r-', linewidth=3.0)
		plt.plot(self.t, self.disp, 'k-')
		plt.xlabel('Time (s)')
		plt.ylabel('Displacement (cm)')
		plt.legend(('Numerical', 'Analytical'))
		plt.savefig('fig_displacement.png')
		plt.close()