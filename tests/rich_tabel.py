import csv
from rich.console import Console
from rich.table import Table

# 初始化一个 Rich Console 对象
console = Console()

# 模拟 CSV 数据
csv_data = [
    ["Name", "Age", "Email"],
    ["Alice", "25", "alice@example.com"],
    ["Bob", "", "bob@example.com"],
    ["", "30", "charlie@example.com"],
    ["David", "35", "david@example.com"],
]

# 初始化成功和失败的计数器
success_count = 0
failure_count = 0

# 模拟文件名
filename = "example.csv"

# 读取 CSV 数据
for row in csv_data[1:]:
    name, age, email = row
    if name and age and email:
        success_count += 1
    else:
        failure_count += 1

# 创建一个表格对象
table = Table(title="CSV Processing Summary")

# 添加表格列
table.add_column("Status", justify="left", style="cyan")
table.add_column("Count", justify="right", style="magenta")
table.add_column("Filename", justify="left", style="green")

# 添加表格行
table.add_row("Success", str(success_count), filename)
table.add_row("Failure", str(failure_count), filename)

# 输出表格
console.print(table)
