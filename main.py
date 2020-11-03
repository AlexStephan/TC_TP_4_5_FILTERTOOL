from src.filterTool import FilterTool
from PyQt5.QtWidgets import QApplication

if __name__ == '__main__':
    app = QApplication([])
    window = FilterTool()
    window.show()
    app.exec()

