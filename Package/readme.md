# Packageの作成方法
- packageとは？
    - moduleのこと、


Template(;
      user="kohei yoshida",
      authors=["kohei yoshida"],
      dir=".",
      julia=v"1.10.4",
      plugins=[
        License(; name="MIT"),
        Git(; ssh=true),
        GitHubActions(),
      ],
  )