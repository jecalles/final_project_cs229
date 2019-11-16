from context import src
from src import test_module

print('passes ../src import')

from context import res
print('passes ../res import')

x = src.test_module.myclass()
x.hello()

