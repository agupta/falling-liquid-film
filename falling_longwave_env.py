import numpy as np
import gym
from gym import spaces
import pexpect

class FlattenLongWaveEnv(gym.Env):
  """
  Custom Environment that follows gym interface.
  This is an env where the agent must learn to give velocities at nozzles such
  that the interface goes flat as quick as possible. Details in writeup that
  isn't written yet :).
  """
  # # Because of google colab, we cannot implement the GUI ('human' render mode)
  # metadata = {'render.modes': ['console']}


  def __init__(self, num_nozzles=8, nozzle_width=7, grid_size=128, scale=0.1):
    super(FlattenLongWaveEnv, self).__init__()

    # Number of nozzles to use for blowing and suction.
    self.num_nozzles = num_nozzles
    self.nozzle_width = nozzle_width
    self.grid_size = grid_size
    self.scale = scale # this scales the velocities.

    # TODO: other params, etc.
    self.solver = None

    # Define action and observation space
    # They must be gym.spaces objects

    # The action will be normalized velocities at each of the first
    # (num_nozzles - 1) nozzles, last one is 0 - the sum of the others to
    # satisfy the condition discussed in ThompsonJFM.
    n_actions = self.num_nozzles - 1
    self.action_space = spaces.Box(low=-1, high=1, shape=(n_actions,), dtype=np.float64)
    # The observation will be the interface height at each grid point.
    self.observation_space = spaces.Box(low=0, high=2, shape=(self.grid_size,), dtype=np.float64)

  def _observe_h(self):
    self.solver.expect_exact("h >")
    index = self.solver.expect_exact(["F >", "finish!"])
    observation = np.fromstring(self.solver.before, dtype=np.float64, sep=' ')
    done = (index == 1) # we observe finish! rather than F >
    return observation, done
  
  def _reward(self, h):
    # simply the negative of the l2 norm of h - 1 (when h is treated like a vector)
    return - np.linalg.norm(h - 1) # NOTE: np.sqrt(x.dot(x)) is faster than np.linalg.norm

  def reset(self):
    """
    Important: the observation must be a numpy array
    :return: (np.array) 
    """
    # # Initialize the agent at the right of the grid
    # self.agent_pos = self.grid_size - 1
    # # here we convert to float32 to make it more general (in case we want to use continuous actions)
    # return np.array([self.agent_pos]).astype(np.float32)
    
    if self.solver is not None:
      self.solver.close()
    self.solver = pexpect.spawn("bin/be_dirty_no_log")
    observation, _ = self._observe_h() # assume can't be done already
    self.last_obs = observation
    return observation

  def _bump(self, x):
    return np.piecewise(x, [(x > 0) & (x < self.nozzle_width), (x >= self.nozzle_width)|(x<=0)], [lambda x: np.exp(1 - 1/(1-((2*x/self.nozzle_width - 1)**2))), 0])
    #return (x>0) * (x<self.nozzle_width) * np.exp(1 - 1/(1-((2*x/self.nozzle_width - 1)**2)))

  def step(self, action):
    action_sum = np.sum(action)

    # Check whether the last nozzle will be out of bounds 
    if not -1 <= action_sum <= 1:
      # don't update solver, very negative reward, return
      observation = self.last_obs
      reward = -1.0 * abs(action_sum)# multiplying makes a gradient so that its easier to learn to make nozzles in bounds.
      done = False
      info = {} # TODO: add info?
      return observation, reward, done, info
    
    # Add the last nozzle to array
    last_nozzle = -action_sum
    nozzle_velocities = np.append(action, last_nozzle)

    # Make the F field and input it
    x = np.arange(128)
    offsets = np.linspace(0, self.grid_size, self.num_nozzles, endpoint=False)
    F = np.zeros(128)
    for nozzle_v, offset in np.column_stack((nozzle_velocities, offsets)):
      F += nozzle_v * self.scale * self._bump(x - offset)
    F_input = np.array2string(F)[1:-1] # remove square brackets
    self.solver.sendline(F_input)

    observation, done = self._observe_h()
    self.last_obs = observation

    reward = 1.0 + self._reward(observation) # add 1 cos it might help who knows

    # Optionally we can pass additional info, we are not using that for now
    info = {}

    return observation, reward, done, info

  def render(self, mode='human'):
    raise NotImplementedError()

  def close(self):
    self.solver.close()
    