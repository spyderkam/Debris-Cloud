## Find average probability of fragments being inside bloated sphere

```python
in_probs_s, out_probs_s = [], []     # inside and outside probabilities of small fragments
in_probs_m, out_probs_m = [], []
in_probs_l, out_probs_l = [], []

for i in range(2):
  smallCloud = Cloud(0.15, 100);  vec_r_s = smallCloud.sample_positions()
  medCloud =  Cloud(0.09, 100);   vec_r_m = medCloud.sample_positions()
  largeCloud = Cloud(0.05, 100);  vec_r_l = largeCloud.sample_positions()
  
  inl_s, outl_s = [], []             # inside and outside numbers of small fragments
  inl_m, outl_m = [], []
  inl_l, outl_l = [], []
  
  for i in range(len(vec_r_s)):
    r_s, r_m, r_l = vec_r_s[i], vec_r_m[i], vec_r_l[i]
    
    fragment_data = {
      'small': (r_s, smallCloud, inl_s, outl_s, in_probs_s, out_probs_s),
      'medium': (r_m, medCloud, inl_m, outl_m, in_probs_m, out_probs_m),
      'large': (r_l, largeCloud, inl_l, outl_l, in_probs_l, out_probs_l),
    }
    
    for size_category in fragment_data.keys():
      r, cloud, inside, outside, inProb, outProb = fragment_data[size_category]
      if np.linalg.norm(r) <= cloud.updated_radius(t=0):
        inside.append(r)
      else:
        outside.append(r)
  
      inProb.append(int(len(inside)/cloud.nFrag*100))
      outProb.append(int(len(outside)/cloud.nFrag*100))

print(f"Small fragments: {np.mean(in_probs_s)}% inside, {np.mean(out_probs_s)}% outside")
print(f"Medium fragments: {np.mean(in_probs_m)}% inside, {np.mean(out_probs_m)}% outside")
print(f"Large fragments: {np.mean(in_probs_l)}% inside, {np.mean(out_probs_l)}% outside")
```


## Basic plot

```python
cloud = Cloud(characteristic_length=0.15, num_fragments=100, breakup_type="collision")
vec_r = init_positions = cloud.sample_positions()

with open('plotter.py', 'r') as file:
    code = file.read()
    exec(code)

inside, outside = [], []
for r in vec_r:
    if np.linalg.norm(r) <= cloud.updated_radius(t=0):
        inside.append(r)
    else:
        outside.append(r)

print(f"{round(len(inside)/cloud.nFrag*100, 2)}% of the fragments are INside.")
print(f"{round(len(outside)/cloud.nFrag*100, 2)}% of the fragments are OUTside.")
```