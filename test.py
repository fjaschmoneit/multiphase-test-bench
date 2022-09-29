

def test(a):
    if a:
        return 1
    else:
        return 0

for i in ['a', 'b', '0', [1,2,3] ]:
    print(test(i))