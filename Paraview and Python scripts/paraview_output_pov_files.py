from paraview.simple import *

timesteps = GetTimeKeeper().TimestepValues
view = GetAnimationScene()
renderView1 = GetActiveViewOrCreate('RenderView')

for i in range(len(timesteps)):
	view.AnimationTime = timesteps[i]
	ExportView('F:/Chyu/rising gate/0408/50/' +str(i+41)+'.pov', view=renderView1)
