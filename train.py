from falling_longwave_env import FlattenLongWaveEnv

from stable_baselines.common.env_checker import check_env
env = FlattenLongWaveEnv()
# If the environment don't follow the interface, an error will be thrown
check_env(env, warn=True)

from stable_baselines import SAC, TD3 # State of the art continuous algorithms

env = FlattenLongWaveEnv()

model = SAC('MlpPolicy', env, verbose=1)
model.learn(total_timesteps=1000, log_interval=10)