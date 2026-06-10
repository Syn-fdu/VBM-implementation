# Version10.1 folder-no-extensions patch

本次修改采用新的结构方案：

- 保留外层 `version10.1/` 文件夹。
- 删除项目根目录下原来的 `extensions/` 文件夹。
- 四个 extension 入口改为根目录下四个独立 R 脚本。
- `main.R` 不再调用 `extensions/run_all_extensions.R`，而是直接 source 四个根目录 extension 脚本。
- enhancedplot Python 脚本统一移动到 `scripts/`，并仍由对应 R 脚本自动调用。
- 输出移动到 `output/version10_extension_results/`。
